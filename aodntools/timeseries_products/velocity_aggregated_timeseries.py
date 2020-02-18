import os
import sys
import tempfile
import shutil
import netCDF4 as nc4
import numpy as np
import json
from datetime import datetime
import argparse
from pkg_resources import resource_filename
from aodntools import __version__

import xarray as xr

import aggregated_timeseries as TStools

TEMPLATE_JSON = resource_filename(__name__, 'velocity_aggregated_timeseries_template.json')


def sort_files(file_list, input_dir=""):
    """
    sort list of files according to deployment date
    :param file_list: List of files to sort
    :return: sorted list of files
    """

    time_start = []
    for file in file_list:
        with nc4.Dataset(os.path.join(input_dir, file), 'r') as ds:
            time_start.append(np.datetime64(ds.time_deployment_start))
    tuples = sorted(zip(time_start, files_to_agg))
    return [t[1] for t in tuples]


def check_file(nc, site_code):
    """
    Return list of errors found in the file:
    Variables of interest are present
    TIME, DEPTH, LATITUDE, LONGITUDE,  is present
    NOMINAL_DEPTH is not present as variable or attribute
    file_version is not FV01
    if LATITUDE or LONIGITUDE dimension has length >1

    :param nc: xarray dataset
    :param site_code: code of the mooring site
    :return: dictionary with the file name and list of failed tests
    """

    attributes = list(nc.attrs)
    variables = list(nc.variables)
    allowed_dimensions = ['TIME', 'LATITUDE', 'LONGITUDE', 'HEIGHT_ABOVE_SENSOR']
    required_variables = ['UCUR', 'VCUR', 'WCUR']
    error_list = []

    if nc.site_code != site_code:
        error_list.append('Wrong site_code: ' + nc.site_code)

    nc_file_version = nc.file_version
    if 'Level 1' not in nc_file_version:
        error_list.append('Wrong file version: ' + nc_file_version)

    if 'DEPTH' not in variables:
        error_list.append('DEPTH variable missing')

    if 'DEPTH' not in variables and 'HEIGHT_ABOVE_SENSOR' not in variables:
        error_list.append('DEPTH and HEIGHT_ABOVE_SENSOR missing')

    if 'TIME' not in variables:
        error_list.append('TIME variable missing')

    if 'LATITUDE' not in variables:
        error_list.append('LATITUDE variable missing')

    if 'LONGITUDE' not in variables:
        error_list.append('LONGITUDE variable missing')


    for variable in required_variables:
        if variable not in variables:
            error_list.append(variable + ' variable missing')
        else:
            VoIdimensions = list(nc[variable].dims)
            if 'TIME' not in VoIdimensions:
                error_list.append('TIME is not a dimension for ' + variable)
            if 'LATITUDE' in VoIdimensions and len(nc.LATITUDE) > 1:
                error_list.append('more than one LATITUDE for ' + variable)
            if 'LONGITUDE' in VoIdimensions and len(nc.LONGITUDE) > 1:
                error_list.append('more than one LONGITUDE for ' + variable)
            for dim in VoIdimensions:
                if dim not in allowed_dimensions:
                    error_list.append('not allowed dimension: ' + dim)

    if 'NOMINAL_DEPTH' not in variables and 'instrument_nominal_depth' not in attributes:
        error_list.append('no NOMINAL_DEPTH')

    return error_list


def get_number_flatvalues(nc):
    """
    Get the number of flatten values and the number of cells above the sensor
    :param nc: xarray dataset
    :return: number of values, number of cells above the sensor
    """
    if 'HEIGHT_ABOVE_SENSOR' in nc.dims:
        n_cells = nc.dims['HEIGHT_ABOVE_SENSOR']
        n_flatt_values = nc.dims['TIME'] * n_cells
    else:
        n_cells = 1
        n_flatt_values = nc.dims['TIME']
    return n_flatt_values, n_cells


def flat_variable(nc, varname):
    """
    Return a 1D array of 2D values
    :param nc: dataset
    :param varname: Variable of interest
    :return: variable values flattened
    """
    return nc[varname].values.flatten()


def get_instrumentID(nc):
    """
    Create instrument id based on deployment metadata
    :param nc: xarray dataset
    :return: instrumentID as string
    """
    return '; '.join([nc.deployment_code, nc.instrument, nc.instrument_serial_number])


def in_water(nc):
    """
    cut data the entire dataset to in-water only timestamps, dropping the out-of-water records.
    :param nc: xarray dataset
    :return: xarray dataset
    """
    time_deployment_start = np.datetime64(nc.attrs['time_deployment_start'][:-1])
    time_deployment_end = np.datetime64(nc.attrs['time_deployment_end'][:-1])
    TIME = nc['TIME'][:]
    return nc.where((TIME >= time_deployment_start) & (TIME <= time_deployment_end), drop=True)




## MAIN FUNCTION
def velocity_aggregated(files_to_agg, site_code, input_dir='', output_dir='./',
                        download_url_prefix=None, opendap_url_prefix=None):
    """
    Aggregate U, V and W CUR variables from all deployments at one site.
    the vertical cells are flattened and related to its depth
    additional metadata variables are stored to track the origin of the data
    :param files_to_agg: list of files to aggregate
    :param site_code: site code
    :param input_dir: base path where source files are stored
    :param output_dir: path where the result file will be written
    :param download_url_prefix: URL prefix for file download (to be prepended to paths in files_to_agg)
    :param opendap_url_prefix: URL prefix for OPENAP access (to be prepended to paths in files_to_agg)
    :return: file path of the aggregated product, list of rejected files
    """

    varlist = ['UCUR', 'VCUR', 'WCUR', 'DEPTH']
    time_units="days since 1950-01-01 00:00:00 UTC"
    time_calendar="gregorian"
    epoch = np.datetime64("1950-01-01T00:00:00")
    one_day = np.timedelta64(1, 'D')

    varlen_list = []
    bad_files = []
    varlen_file = []
    rejected_files = []

    # default name for temporary file. It will be renamed at the end
    _, temp_outfile = tempfile.mkstemp(suffix='.nc', dir=output_dir)

    ## sort the file list in chronological order
    files_to_agg = sort_files(files_to_agg, input_dir=input_dir)

    ## check files and get total number of flattened obs
    for file in files_to_agg:
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            ## clip to in water data only
            nc = in_water(nc)

            varlen_file.append(get_number_flatvalues(nc))
            error_list = check_file(nc, site_code)
            if not error_list:
                varlen_list.append(get_number_flatvalues(nc)[0])
            else:
                bad_files.append([file, error_list])
                rejected_files.append(file)

    ## remove bad files form the list
    for file in bad_files:
        files_to_agg.remove(file[0])


    varlen_list = [0] + varlen_list
    varlen_total = sum(varlen_list)
    n_files = len(files_to_agg)


    ## create ncdf file, dimensions and variables
    ds = nc4.Dataset(os.path.join(output_dir, temp_outfile), 'w')
    OBSERVATION = ds.createDimension('OBSERVATION', size=varlen_total)
    INSTRUMENT = ds.createDimension('INSTRUMENT', size=n_files)

    obs_double_template = {'datatype': np.float64, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_float_template = {'datatype': np.float32, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_byte_template = {'datatype': 'byte', 'zlib': True, 'dimensions': ('OBSERVATION'), 'fill_value': 99}
    obs_int_template = {'datatype': 'int', 'zlib': True, 'dimensions': ('OBSERVATION')}
    inst_S256_template = {'datatype': 'str', 'dimensions': ('INSTRUMENT')}
    inst_float_template ={'datatype': np.float32, 'dimensions': ('INSTRUMENT')}
    inst_double_template ={'datatype': np.float64, 'dimensions': ('INSTRUMENT')}



    UCUR = ds.createVariable(varname='UCUR', **obs_float_template)
    VCUR = ds.createVariable(varname='VCUR', **obs_float_template)
    WCUR = ds.createVariable(varname='WCUR', **obs_float_template)
    DEPTH = ds.createVariable(varname='DEPTH', **obs_float_template)
    UCURqc = ds.createVariable(varname='UCUR_quality_control', **obs_byte_template)
    VCURqc = ds.createVariable(varname='VCUR_quality_control', **obs_byte_template)
    WCURqc = ds.createVariable(varname='WCUR_quality_control', **obs_byte_template)
    DEPTHqc = ds.createVariable(varname='DEPTH_quality_control', **obs_byte_template)
    TIME = ds.createVariable(varname='TIME', **obs_double_template)
    instrument_index = ds.createVariable(varname='instrument_index', **obs_int_template)

    source_file = ds.createVariable(varname='source_file', **inst_S256_template)
    instrument_id = ds.createVariable(varname='instrument_id', **inst_S256_template)
    LATITUDE = ds.createVariable(varname='LATITUDE', **inst_double_template)
    LONGITUDE = ds.createVariable(varname='LONGITUDE', **inst_double_template)
    NOMINAL_DEPTH = ds.createVariable(varname='NOMINAL_DEPTH', **inst_float_template)
    SECONDS_TO_MIDDLE = ds.createVariable(varname='SECONDS_TO_MIDDLE', **inst_float_template)

    ## main loop
    for index, file in enumerate(files_to_agg):
        print(index)
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:

            ## in water data only
            nc = in_water(nc)

            start = sum(varlen_list[:index + 1])
            end = sum(varlen_list[:index + 2])
            n_cells = get_number_flatvalues(nc)[1]
            UCUR[start:end] = flat_variable(nc, 'UCUR')
            UCURqc[start:end] = flat_variable(nc, 'UCUR_quality_control')
            VCUR[start:end] = flat_variable(nc, 'VCUR')
            VCURqc[start:end] = flat_variable(nc, 'VCUR_quality_control')
            if 'WCUR' in nc.data_vars:
                WCUR[start:end] = flat_variable(nc, 'WCUR')
                WCURqc[start:end] = flat_variable(nc, 'WCUR_quality_control')
            else:
                WCUR[start:end] = np.full(varlen_list[index], np.nan)
                WCURqc[start:end] = np.full(varlen_list[index], np.nan)
            ##calculate depth
            if 'HEIGHT_ABOVE_SENSOR' in nc.dims:
                DEPTH[start:end] = (nc.DEPTH - nc.HEIGHT_ABOVE_SENSOR).values.flatten()
                DEPTHqc[start:end] = np.array(n_cells * [nc.DEPTH_quality_control.values]).flatten()
            else:
                DEPTH[start:end] = nc.DEPTH.values
                DEPTHqc[start:end] = nc.DEPTH_quality_control.values
            ## set TIME and instrument index
            TIME[start:end] = (np.repeat(flat_variable(nc, 'TIME'), n_cells) - epoch) / one_day
            instrument_index[start:end] = np.repeat(index, varlen_list[index + 1])
            ## get and store deployment metadata
            LATITUDE[index] = nc.LATITUDE.values
            LONGITUDE[index] = nc.LONGITUDE.values
            NOMINAL_DEPTH[index] = TStools.get_nominal_depth(nc)
            source_file[index] = file
            instrument_id[index] = get_instrumentID(nc)
            ## add time offset to the middle of the measuring window, if it exists
            if 'seconds_to_middle_of_measurement' in nc.TIME.attrs:
                SECONDS_TO_MIDDLE[index] = nc.TIME.seconds_to_middle_of_measurement
            else:
                SECONDS_TO_MIDDLE[index] = np.nan

    ## add atributes
    with open(TEMPLATE_JSON) as json_file:
        attribute_dictionary = json.load(json_file)
    variable_attribute_dictionary = attribute_dictionary['_variables']
    global_attribute_dictionary = attribute_dictionary['_global']

    ## set variable attrs
    for var in list(ds.variables):
        ds[var].setncatts(variable_attribute_dictionary[var])

    if download_url_prefix or opendap_url_prefix:
        ds['source_file'].setncatts(TStools.source_file_attributes(download_url_prefix, opendap_url_prefix))

    ## set global attrs
    timeformat = '%Y-%m-%dT%H:%M:%SZ'
    file_timeformat = '%Y%m%d'

    time_start = nc4.num2date(np.min(TIME[:]), time_units, time_calendar).strftime(timeformat)
    time_end = nc4.num2date(np.max(TIME[:]), time_units, time_calendar).strftime(timeformat)
    time_start_filename = nc4.num2date(np.min(TIME[:]), time_units, time_calendar).strftime(file_timeformat)
    time_end_filename = nc4.num2date(np.max(TIME[:]), time_units, time_calendar).strftime(file_timeformat)

    contributor_name, contributor_email, contributor_role = TStools.get_contributors(files_to_agg=files_to_agg, input_dir=input_dir)
    add_attribute = {
                    'title':                    ("Long Timeseries Velocity Aggregated product: " + ', '.join(varlist) + " at " +
                                                  site_code + " between " + time_start + " and " + time_end),
                    'site_code':                site_code,
                    'time_coverage_start':      time_start,
                    'time_coverage_end':        time_end,
                    'geospatial_vertical_min':  np.min(ds['DEPTH']),
                    'geospatial_vertical_max':  np.max(ds['DEPTH']),
                    'geospatial_lat_min':       np.min(ds['LATITUDE']),
                    'geospatial_lat_max':       np.max(ds['LATITUDE']),
                    'geospatial_lon_min':       np.min(ds['LONGITUDE']),
                    'geospatial_lon_max':       np.max(ds['LONGITUDE']),
                    'date_created':             datetime.utcnow().strftime(timeformat),
                    'history':                  datetime.utcnow().strftime(timeformat) + ': Aggregated file created.',
                    'keywords':                 ', '.join(varlist + ['AGGREGATED']),
                    'rejected_files':           "\n".join(rejected_files),
                    'contributor_name':        "; ".join(contributor_name),
                    'contributor_email':       "; ".join(contributor_email),
                    'contributor_role':        "; ".join(contributor_role),
                    'generating_code_version':  __version__
    }

    ## add version
    github_comment = ('\nThis file was created using https://github.com/aodn/python-aodntools/blob/'
                      '{v}/aodntools/timeseries_products/aggregated_timeseries.py'.format(v=__version__)
                      )
    global_attribute_dictionary['lineage'] += github_comment

    global_attribute_dictionary.update(add_attribute)
    ds.setncatts(dict(sorted(global_attribute_dictionary.items())))

    ds.close()


    ## create the output file name and rename the tmp file
    facility_code = TStools.get_facility_code(os.path.join(input_dir, files_to_agg[0]))
    data_code = 'VZ'
    product_type = 'aggregated-timeseries'
    file_version = 1
    output_name = '_'.join(['IMOS', facility_code, data_code, time_start_filename, site_code, ('FV0'+str(file_version)),
                            ("velocity-"+product_type),
                            ('END-'+ time_end_filename), 'C-' + datetime.utcnow().strftime(file_timeformat)]) + '.nc'
    ncout_path = os.path.join(output_dir, output_name)
    shutil.move(temp_outfile, ncout_path)


    return ncout_path, bad_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenate X,Y,Z velocity variables from ALL instruments from ALL deployments from ONE site")
    parser.add_argument('-site', dest='site_code', help='site code, like NRMMAI',  required=True)
    parser.add_argument('-files', dest='filenames', help='name of the file that contains the source URLs', required=True)
    parser.add_argument('-path', dest='output_path', help='path where the result file will be written. Default: ./', default='./', required=False)
    parser.add_argument('-indir', dest='input_dir', help='base path of input files', default='', required=False)
    parser.add_argument('-outdir', dest='output_dir', help='path where the result file will be written. Default ./',
                        default='./', required=False)
    parser.add_argument('-download_url', dest='download_url', help='path to the download_url_prefix',
                        default='', required=False)
    parser.add_argument('-opendap_url', dest='opendap_url', help='path to the opendap_url_prefix',
                        default='', required=False)

    args = parser.parse_args()

    with open(args.filenames) as ff:
        files_to_agg = [line.rstrip() for line in ff]

    print(velocity_aggregated(files_to_agg=files_to_agg, site_code=args.site_code,
                              input_dir=args.input_dir, output_dir=args.output_dir,
                              download_url_prefix=args.download_url, opendap_url_prefix=args.opendap_url))
