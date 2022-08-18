#!/usr/bin/env python3

import argparse
from copy import deepcopy
import json
import os
import shutil
import tempfile

import numpy as np
import xarray as xr
from netCDF4 import Dataset, num2date, stringtochar
from pkg_resources import resource_filename

from aodntools import __version__
from aodntools.timeseries_products.common import (
    NoInputFilesError,
    check_file,
    in_water,
    current_utc_timestamp,
    TIMESTAMP_FORMAT,
    DATESTAMP_FORMAT,
)

TEMPLATE_JSON = resource_filename(__name__, "aggregated_timeseries_template.json")


def sort_files(files_to_agg, input_dir=""):
    """
    sort list of files according to deployment date
    :param files_to_agg: List of files to sort
    :param input_dir: base path where source files are stored
    :return: sorted list of files
    """

    time_start = []
    for file in files_to_agg:
        with Dataset(os.path.join(input_dir, file)) as ds:
            time_start.append(np.datetime64(ds.time_deployment_start.rstrip("Z")))
    tuples = sorted(zip(time_start, files_to_agg))
    return [t[1] for t in tuples]


def get_variable_values(nc, variable):
    """
    Get values of the variable and its QC flags.
    If variable is not present, nan returned, its QC flags set to 9
    If variable present but not its QC flags, QC set to 0
    :param nc: xarray dataset
    :param variable: name of the variable to get
    :return: variable values and variable qc flags
    """
    n_records = len(nc.TIME)
    file_variables = list(nc.variables)

    if variable in file_variables:
        variable_values = nc[variable].values
        if any(np.isnan(variable_values)):
            variable_values = np.ma.masked_array(
                variable_values, mask=np.isnan(variable_values)
            )
        if variable + "_quality_control" in file_variables:
            variableQC_values = nc[variable + "_quality_control"].values
        else:
            variableQC_values = 0
    else:
        variable_values = np.ma.masked
        variableQC_values = 9

    return variable_values, variableQC_values


def get_instrument_id(nc):
    """
    Create instrument id based on deployment metadata
    :param nc: xarray dataset
    :return: instrumentID as string
    """
    deployment_code = nc.attrs.get("deployment_code", "[unknown deployment]")
    instrument = nc.attrs.get("instrument", "[unknown instrument]")
    instrument_serial_number = nc.attrs.get(
        "instrument_serial_number", "[unknown serial number]"
    )
    return "; ".join([deployment_code, instrument, instrument_serial_number])


def get_nominal_depth(nc):
    """
    retunr nominal depth from NOMINAL_DEPTH variable or
    if it is not present from instrument_nominal_depth global attribute

    :param nc: xarray dataset
    :return: nominal depth of the instrument
    """

    if "NOMINAL_DEPTH" in list(nc.variables):
        nominal_depth = nc.NOMINAL_DEPTH.squeeze().values
    else:
        nominal_depth = nc.instrument_nominal_depth

    return nominal_depth


def get_contributors(files_to_agg, input_dir=""):
    """
    get the author and principal investigator details for each file

    :param files_to_agg: list of files
    :param input_dir: base path where source files are stored
    :return: dict of attributes contributor_name, contributor_email and contributor_role
    """

    contributors = set()
    contributor_name, contributor_email, contributor_role = [], [], []

    for file in files_to_agg:
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            attributes = nc.attrs.keys()
            if all(att in attributes for att in ["author", "author_email"]):
                contributors.add((nc.author, nc.author_email, "author"))
            if all(
                att in attributes
                for att in ["principal_investigator", "principal_investigator_email"]
            ):
                contributors.add(
                    (
                        nc.principal_investigator,
                        nc.principal_investigator_email,
                        "principal_investigator",
                    )
                )
    for item in contributors:
        contributor_name.append(item[0])
        contributor_email.append(item[1])
        contributor_role.append(item[2])

    return {
        "contributor_name": "; ".join(contributor_name),
        "contributor_email": "; ".join(contributor_email),
        "contributor_role": "; ".join(contributor_role),
    }


def get_data_code(VoI):
    """
    get data code sensu IMOS conventions from variable code

    :param VoI: variable code
    :return: variable data code
    """

    # dictionary of data code. could be read from external file
    dataCodes = {
        "DEPTH": "Z",
        "PRES": "Z",
        "PRES_REL": "Z",
        "TEMP": "T",
        "PSAL": "S",
        "PAR": "F",
        "TURB": "U",
        "TURBF": "U",
        "DOX1": "O",
        "DOX1_2": "O",
        "DOX1_3": "O",
        "DOX2": "O",
        "DOX2_1": "O",
        "DOXS": "O",
        "CPHL": "B",
        "CHLU": "B",
        "CHLF": "B",
        "UCUR": "V",
        "VCUR": "V",
        "WCUR": "V",
    }
    return dataCodes[VoI]


def get_facility_code(fileURL):
    """
    get the facility code from the file URL

    :param fileURL: URL of a file
    :return: facility code
    """

    return os.path.basename(fileURL).split("_")[1]


def source_file_attributes(download_url_prefix, opendap_url_prefix):
    """
    If relevant URL prefixes are specified, return attributes to add to the source_file variable to describe their use.
    :param download_url_prefix: prefix string for download URLs
    :param opendap_url_prefix: prefix string for OPENDAP URLs
    :return: dictionary of attributes to add to the source_file variable
    """
    attributes = {
        "comment": "This variable lists the relative path of each input file."
    }
    if download_url_prefix:
        attributes["comment"] += (
            " To obtain a download URL for a file, "
            "append its path to the download_url_prefix attribute."
        )
        attributes["download_url_prefix"] = download_url_prefix
    if opendap_url_prefix:
        attributes["comment"] += (
            " To interact with the file remotely via the OPENDAP protocol, "
            "append its path to the opendap_url_prefix attribute."
        )
        attributes["opendap_url_prefix"] = opendap_url_prefix
    return attributes


## MAIN FUNCTION
def main_aggregator(
    files_to_agg,
    var_to_agg,
    site_code,
    input_dir="",
    output_dir="./",
    download_url_prefix=None,
    opendap_url_prefix=None,
):
    """
    Aggregate the Variable of Interest (VoI) from all deployments at one site.
    additional metadata variables are stored to track the origin of the data
    :param files_to_agg: List of files to aggregate. Each path is interpreted relative to input_dir (if specified).
                         These relative paths are listed in the `source_files` variable in the output file.
    :param site_code: site code
    :param var_to_agg: Variable of Interest
    :param input_dir: base path where source files are stored
    :param output_dir: path where the result file will be written
    :param download_url_prefix: URL prefix for file download (to be prepended to paths in files_to_agg)
    :param opendap_url_prefix: URL prefix for OPENAP access (to be prepended to paths in files_to_agg)
    :return: name of the resulting file, list of rejected files
    """

    time_units = "days since 1950-01-01 00:00:00 UTC"
    time_calendar = "gregorian"
    epoch = np.datetime64("1950-01-01T00:00:00")
    one_day = np.timedelta64(1, "D")

    bad_files = {}
    rejected_files = []

    # default name for temporary file. It will be renamed at the end
    _, temp_outfile = tempfile.mkstemp(suffix=".nc", dir=output_dir)

    ## check files and get total number of flattened obs
    n_obs_total = 0
    for file in files_to_agg:
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:

            error_list = check_file(nc, site_code, var_to_agg)
            if not error_list:
                nc = in_water(nc)
                n_obs_total += len(nc.TIME)
            else:
                bad_files.update({file: error_list})
                rejected_files.append(file)

    ## remove bad files form the list and sort in chronological order
    files_to_agg = deepcopy(files_to_agg)
    for file in bad_files.keys():
        files_to_agg.remove(file)
    if len(files_to_agg) == 0:
        raise NoInputFilesError("no valid input files to aggregate")
    files_to_agg = sort_files(files_to_agg, input_dir=input_dir)

    n_files = len(files_to_agg)

    ## create ncdf file, dimensions and variables
    ds = Dataset(os.path.join(output_dir, temp_outfile), "w", format="NETCDF4_CLASSIC")
    OBSERVATION = ds.createDimension("OBSERVATION", size=n_obs_total)
    INSTRUMENT = ds.createDimension("INSTRUMENT", size=n_files)
    STRING256 = ds.createDimension("strlen", 256)

    obs_float_template = {
        "datatype": np.float32,
        "zlib": True,
        "dimensions": ("OBSERVATION",),
        "fill_value": 99999.0,
    }
    obs_double_template = {
        "datatype": np.float64,
        "zlib": True,
        "dimensions": ("OBSERVATION",),
        "fill_value": 99999.0,
    }
    obs_byte_template = {
        "datatype": np.byte,
        "zlib": True,
        "dimensions": ("OBSERVATION",),
        "fill_value": 99,
    }
    obs_int_template = {
        "datatype": np.int16,
        "zlib": True,
        "dimensions": ("OBSERVATION",),
    }
    inst_S256_template = {"datatype": "S1", "dimensions": ("INSTRUMENT", "strlen")}
    inst_float_template = {
        "datatype": np.float32,
        "dimensions": ("INSTRUMENT",),
        "fill_value": 99999.0,
    }
    inst_double_template = {
        "datatype": np.float64,
        "dimensions": ("INSTRUMENT",),
        "fill_value": 99999.0,
    }

    agg_variable = ds.createVariable(varname=var_to_agg, **obs_float_template)
    agg_variable_qc = ds.createVariable(
        varname=var_to_agg + "_quality_control", **obs_byte_template
    )
    DEPTH = ds.createVariable(varname="DEPTH", **obs_float_template)
    DEPTHqc = ds.createVariable(varname="DEPTH_quality_control", **obs_byte_template)
    PRES = ds.createVariable(varname="PRES", **obs_float_template)
    PRESqc = ds.createVariable(varname="PRES_quality_control", **obs_byte_template)
    PRES_REL = ds.createVariable(varname="PRES_REL", **obs_float_template)
    PRES_RELqc = ds.createVariable(
        varname="PRES_REL_quality_control", **obs_byte_template
    )

    TIME = ds.createVariable(varname="TIME", **obs_double_template)
    instrument_index = ds.createVariable(varname="instrument_index", **obs_int_template)

    source_file = ds.createVariable(varname="source_file", **inst_S256_template)
    instrument_id = ds.createVariable(varname="instrument_id", **inst_S256_template)
    LATITUDE = ds.createVariable(varname="LATITUDE", **inst_double_template)
    LONGITUDE = ds.createVariable(varname="LONGITUDE", **inst_double_template)
    NOMINAL_DEPTH = ds.createVariable(varname="NOMINAL_DEPTH", **inst_float_template)

    ## main loop
    start = 0
    for index, file in enumerate(files_to_agg):
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            nc = in_water(nc)
            n_obs = len(nc.TIME)
            end = start + n_obs
            agg_variable[start:end], agg_variable_qc[start:end] = get_variable_values(
                nc, var_to_agg
            )
            DEPTH[start:end], DEPTHqc[start:end] = get_variable_values(nc, "DEPTH")
            PRES[start:end], PRESqc[start:end] = get_variable_values(nc, "PRES")
            PRES_REL[start:end], PRES_RELqc[start:end] = get_variable_values(
                nc, "PRES_REL"
            )

            ## set TIME and instrument index
            TIME[start:end] = (nc.TIME.values - epoch) / one_day
            instrument_index[start:end] = np.repeat(index, n_obs)
            ## get and store deployment metadata
            LATITUDE[index] = nc.LATITUDE.values
            LONGITUDE[index] = nc.LONGITUDE.values
            NOMINAL_DEPTH[index] = get_nominal_depth(nc)
            source_file[index] = stringtochar(np.array(file, dtype="S256"))
            instrument_id[index] = stringtochar(
                np.array(get_instrument_id(nc), dtype="S256")
            )

        start = end

    ## add atributes
    with open(TEMPLATE_JSON) as json_file:
        attribute_dictionary = json.load(json_file)
    variable_attribute_dictionary = attribute_dictionary["_variables"]
    global_attribute_dictionary = attribute_dictionary["_global"]

    ## set variable attrs
    for var in list(ds.variables):
        ds[var].setncatts(variable_attribute_dictionary[var])

    if download_url_prefix or opendap_url_prefix:
        ds["source_file"].setncatts(
            source_file_attributes(download_url_prefix, opendap_url_prefix)
        )

    ## set global attrs
    time_start = num2date(np.min(TIME[:]), time_units, time_calendar).strftime(
        TIMESTAMP_FORMAT
    )
    time_end = num2date(np.max(TIME[:]), time_units, time_calendar).strftime(
        TIMESTAMP_FORMAT
    )
    time_start_filename = num2date(np.min(TIME[:]), time_units, time_calendar).strftime(
        DATESTAMP_FORMAT
    )
    time_end_filename = num2date(np.max(TIME[:]), time_units, time_calendar).strftime(
        DATESTAMP_FORMAT
    )

    add_attribute = {
        "title": (
            "Long Timeseries Velocity Aggregated product: "
            + var_to_agg
            + " at "
            + site_code
            + " between "
            + time_start
            + " and "
            + time_end
        ),
        "site_code": site_code,
        "time_coverage_start": time_start,
        "time_coverage_end": time_end,
        "geospatial_vertical_min": np.min(ds["DEPTH"][:]),
        "geospatial_vertical_max": np.max(ds["DEPTH"][:]),
        "geospatial_lat_min": np.min(ds["LATITUDE"][:]),
        "geospatial_lat_max": np.max(ds["LATITUDE"][:]),
        "geospatial_lon_min": np.min(ds["LONGITUDE"][:]),
        "geospatial_lon_max": np.max(ds["LONGITUDE"][:]),
        "date_created": current_utc_timestamp(),
        "history": current_utc_timestamp() + ": Aggregated file created.",
        "keywords": ", ".join([var_to_agg, "AGGREGATED"]),
        "rejected_files": "\n".join(rejected_files),
        "generating_code_version": __version__,
    }
    add_attribute.update(
        get_contributors(files_to_agg=files_to_agg, input_dir=input_dir)
    )

    github_comment = (
        "\nThis file was created using https://github.com/aodn/python-aodntools/blob/"
        "{v}/aodntools/timeseries_products/aggregated_timeseries.py".format(
            v=__version__
        )
    )
    global_attribute_dictionary["lineage"] += github_comment
    global_attribute_dictionary.update(add_attribute)
    ds.setncatts(dict(sorted(global_attribute_dictionary.items())))

    ds.close()

    ## create the output file name and rename the tmp file
    facility_code = get_facility_code(os.path.join(input_dir, files_to_agg[0]))
    data_code = get_data_code(var_to_agg) + "Z"
    product_type = "aggregated-timeseries"
    file_version = 1
    output_name = (
        "_".join(
            [
                "IMOS",
                facility_code,
                data_code,
                time_start_filename,
                site_code,
                ("FV0" + str(file_version)),
                (var_to_agg + "-" + product_type),
                ("END-" + time_end_filename),
                "C-" + current_utc_timestamp(DATESTAMP_FORMAT),
            ]
        )
        + ".nc"
    )
    ncout_path = os.path.join(output_dir, output_name)
    shutil.move(temp_outfile, os.path.join(output_dir, ncout_path))

    return ncout_path, bad_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Aggregate ONE variable from ALL instruments from ALL deployments from ONE site"
    )
    parser.add_argument(
        "-site", dest="site_code", help="site code, like NRMMAI", required=True
    )
    parser.add_argument(
        "-var", dest="varname", help="variable to aggregate, like TEMP", required=True
    )
    parser.add_argument(
        "-files",
        dest="filenames",
        help="name of the file that contains the source URLs (relative to inpath, if given)",
        required=True,
    )
    parser.add_argument(
        "-indir",
        dest="input_dir",
        help="base path of input files",
        default="",
        required=False,
    )
    parser.add_argument(
        "-outdir",
        dest="output_dir",
        help="path where the result file will be written. Default ./",
        default="./",
        required=False,
    )
    parser.add_argument(
        "-download_url",
        dest="download_url",
        help="path to the download_url_prefix",
        default="",
        required=False,
    )
    parser.add_argument(
        "-opendap_url",
        dest="opendap_url",
        help="path to the opendap_url_prefix",
        default="",
        required=False,
    )
    args = parser.parse_args()

    with open(os.path.join(args.input_dir, args.filenames)) as ff:
        files_to_agg = [line.rstrip() for line in ff]

    print(
        main_aggregator(
            files_to_agg=files_to_agg,
            var_to_agg=args.varname,
            site_code=args.site_code,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            download_url_prefix=args.download_url,
            opendap_url_prefix=args.opendap_url,
        )
    )
