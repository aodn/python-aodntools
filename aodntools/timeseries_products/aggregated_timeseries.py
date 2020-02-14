#!/usr/bin/env python3

import os
import sys
import uuid
import netCDF4 as nc4
import numpy as np
import json
from datetime import datetime
import argparse
from pkg_resources import resource_filename

import xarray as xr

from aodntools import __version__

TEMPLATE_JSON = resource_filename(__name__, 'aggregated_timeseries_template.json')


def sort_files(files_to_agg):
    """
    sort list of files according to deployment date
    requires netcdf4 as nc4, dateutil.parser as parse
    :param files_to_agg: List of files to sort
    :return: sorted list of files
    """

    time_start = []
    for file in files_to_agg:
        #print(file)
        print('.', end='', flush=True)
        with nc4.Dataset(file, 'r') as ds:
            time_start.append(np.datetime64(ds.time_deployment_start))
    tuples = sorted(zip(time_start, files_to_agg))
    return [t[1] for t in tuples]


def check_file(nc, site_code, VoI):
    """
    Return list of errors found in the file:
    Variables of interest are present
    TIME. LATITUDE, LONGITUDE,  is present
    NOMINAL_DEPTH is not present as variable or attribute
    file_version is not FV01
    if LATITUDE or LONIGUTDE dimension has length >1

    :param nc: xarray dataset
    :param site_code: code of the mooring site
    :param VoI: string. Variable of Interest
    :return: dictionary with the file name and list of failed tests
    """

    attributes = list(nc.attrs)
    file_variables = list(nc.variables)
    allowed_dimensions = ['TIME', 'LATITUDE', 'LONGITUDE']
    error_list = []

    nc_site_code = nc.site_code
    if nc_site_code != site_code:
        error_list.append('Wrong site_code: ' + nc_site_code)

    nc_file_version = nc.file_version
    if 'Level 1' not in nc_file_version:
        error_list.append('Wrong file version: ' + nc_file_version)

    if VoI not in file_variables:
        error_list.append(VoI + ' not in file')
    else:
        for dimension in list(nc[VoI].dims):
            if dimension not in allowed_dimensions:
                error_list.append(dimension+' is not an allowed dimension for ' + VoI)

    if 'NOMINAL_DEPTH' not in file_variables and 'instrument_nominal_depth' not in attributes:
        error_list.append('no NOMINAL_DEPTH')

    return error_list


def get_variable_values(nc, variable):
    """
    Get value sof the variable and its QC flags.
    If variable is not present, nan returned
    If variable present but not its QC flags, QC set to 9
    :param nc: dataset
    :param variable: name of the variable to get
    :return: variable values and variable qc flags
    """
    n_records = len(nc.TIME)
    file_variables = list(nc.variables)

    if variable in file_variables:
        variable_values = nc[variable]
        if variable+'_quality_control' in file_variables:
            variableQC_values = nc[variable+'_quality_control']
        else:
            variableQC_values = np.repeat(9, n_records)
    else:
        variable_values = np.repeat(np.nan, n_records)
        variableQC_values = np.repeat(np.nan, n_records)

    return variable_values, variableQC_values


def get_varvalues(nc, varname):
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


def get_nominal_depth(nc):
    """
    retunr nominal depth from NOMINAL_DEPTH variable or
    if it is not present from instrument_nominal_depth global attribute

    :param nc: xarray dataset
    :return: nominal depth of the instrument
    """

    if 'NOMINAL_DEPTH' in list(nc.variables):
        nominal_depth = nc.NOMINAL_DEPTH.squeeze().values
    else:
        nominal_depth = nc.instrument_nominal_depth

    return nominal_depth


def get_contributors(files_to_agg, input_dir=''):
    """
    get the author and principal investigator details for each file

    :param files_to_aggregate: list of files
    :param input_dir: base path where source files are stored
    :return: list: contributor_name, email and role
    """

    contributors = set()
    contributor_name, contributor_email, contributor_role = [], [], []

    for file in files_to_agg:
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            attributes = nc.attrs.keys()
            if all(att in attributes for att in ['author', 'author_email']):
                contributors.add((nc.author, nc.author_email, 'author'))
            if all(att in attributes for att in ['principal_investigator', 'principal_investigator_email']):
                contributors.add((nc.principal_investigator, nc.principal_investigator_email, 'principal_investigator'))
    for item in contributors:
        contributor_name.append(item[0])
        contributor_email.append(item[1])
        contributor_role.append(item[2])

    return contributor_name, contributor_email, contributor_role


def get_data_code(VoI):
    """
    get data code sensu IMOS conventions from variable code

    :param VoI: variable code
    :return: variable data code
    """

    #dictionary of data code. could be read from external file
    dataCodes = {'DEPTH':       'Z',
                 'PRES':        'Z',
                 'PRES_REL':    'Z',
                 'TEMP':        'T',
                 'PSAL':        'S',
                 'PAR':         'F',
                 'TURB':        'U',
                 'TURBF':       'U',
                 'DOX1':        'O',
                 'DOX1_2':      'O',
                 'DOX1_3':      'O',
                 'DOX2':        'O',
                 'DOX2_1':      'O',
                 'DOXS':        'O',
                 'CPHL':        'B',
                 'CHLU':        'B',
                 'CHLF':        'B',
                 'UCUR':        'V',
                 'VCUR':        'V',
                 'WCUR':        'V'}
    return dataCodes[VoI]


def get_facility_code(fileURL):
    """
    get the facility code from the file URL

    :param fileURL: URL of a file
    :return: facility code
    """

    return os.path.basename(fileURL).split("_")[1]



## MAIN FUNCTION
def aggregate_timeseries(files_to_agg, site_code, VoI, base_path):
    """
    Aggregate the Variable of Interest (VoI) from all deployments at one site.
    additional metadata variables are stored to track the origin of the data
    :param files_to_agg: list of files to aggregate
    :param site_code: site code
    :param VoI: Variable of Interest
    :param base_path: base path to store the resulting file
    :return: name of the resulting file, list of rejected files
    """

    time_units="days since 1950-01-01 00:00:00 UTC"
    time_calendar="gregorian"
    epoch = np.datetime64("1950-01-01T00:00:00")
    one_day = np.timedelta64(1, 'D')

    varlen_list = []
    bad_files = []
    rejected_files = []

    # default name for temporary file. It will be renamed at the end
    outfile = str(uuid.uuid4().hex) + '.nc'

    ## sort the file list in chronological order
    print("SORTING FILES")
    files_to_agg = sort_files(files_to_agg)

    ## check files and get total number of flattened obs
    print('CHECK FILES')
    for file in files_to_agg:
        print('.', end="", flush=True)
        with xr.open_dataset(file) as nc:
            ## clip to in water data only
            nc = in_water(nc)

            error_list = check_file(nc, site_code, VoI)
            if not error_list:
                varlen_list.append(len(nc.TIME))
            else:
                bad_files.append([file, error_list])
                rejected_files.append(file)

    #print(bad_files)

    ## remove bad files form the list
    for file in bad_files:
        files_to_agg.remove(file[0])


    varlen_list = [0] + varlen_list
    varlen_total = sum(varlen_list)
    n_files = len(files_to_agg)

    ## create ncdf file, dimensions and variables
    ds = nc4.Dataset(os.path.join(base_path, outfile), 'w')
    OBSERVATION = ds.createDimension('OBSERVATION', size=varlen_total)
    INSTRUMENT = ds.createDimension('INSTRUMENT', size=n_files)

    obs_float_template = {'datatype': 'float', 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_byte_template = {'datatype': 'byte', 'zlib': True, 'dimensions': ('OBSERVATION'), 'fill_value': 99}
    obs_int_template = {'datatype': 'int', 'zlib': True, 'dimensions': ('OBSERVATION')}
    inst_S256_template = {'datatype': 'S256', 'dimensions': ('INSTRUMENT')}
    inst_float_template ={'datatype': 'float', 'dimensions': ('INSTRUMENT'), "fill_value": 99999.0}



    agg_variable = ds.createVariable(varname=VoI, **obs_float_template)
    agg_variable_qc = ds.createVariable(varname=VoI+'_quality_control', **obs_byte_template)
    DEPTH = ds.createVariable(varname='DEPTH', **obs_float_template)
    DEPTHqc = ds.createVariable(varname='DEPTH_quality_control', **obs_byte_template)
    PRES = ds.createVariable(varname='PRES', **obs_float_template)
    PRESqc = ds.createVariable(varname='PRES_quality_control', **obs_byte_template)
    PRES_REL = ds.createVariable(varname='PRES_REL', **obs_float_template)
    PRES_RELqc = ds.createVariable(varname='PRES_REL_quality_control', **obs_byte_template)

    TIME = ds.createVariable(varname='TIME', **obs_float_template)
    instrument_index = ds.createVariable(varname='instrument_index', **obs_int_template)

    source_file = ds.createVariable(varname='source_file', **inst_S256_template)
    instrument_id = ds.createVariable(varname='instrument_id', **inst_S256_template)
    LATITUDE = ds.createVariable(varname='LATITUDE', **inst_float_template)
    LONGITUDE = ds.createVariable(varname='LONGITUDE', **inst_float_template)
    NOMINAL_DEPTH = ds.createVariable(varname='NOMINAL_DEPTH', **inst_float_template)

    ## main loop
    print('AGGREGATE FILES')
    for index, file in enumerate(files_to_agg):
        print(index, end=",", flush=True)
        with xr.open_dataset(file) as nc:
            nc = in_water(nc)
            file_variables = list(nc.variables)
            start = sum(varlen_list[:index + 1])
            end = sum(varlen_list[:index + 2])
            agg_variable[start:end], agg_variable_qc[start:end] = get_variable_values(nc, VoI)
            DEPTH[start:end], DEPTHqc[start:end] = get_variable_values(nc, 'DEPTH')
            PRES[start:end], PRESqc[start:end] = get_variable_values(nc, 'PRESS')
            PRES_REL[start:end], PRES_RELqc[start:end] = get_variable_values(nc, 'PRESS_REL')

            ## set TIME and instrument index
            TIME[start:end] = (nc.TIME.values - epoch) / one_day
            ##TIME[start:end] = nc.TIME
            instrument_index[start:end] = np.repeat(index, varlen_list[index+1])
            ## get and store deployment metadata
            LATITUDE[index] = nc.LATITUDE.values
            LONGITUDE[index] = nc.LONGITUDE.values
            NOMINAL_DEPTH[index] = get_nominal_depth(nc)
            source_file[index] = file
            instrument_id[index] = get_instrumentID(nc)


    ## add atributes
    with open(TEMPLATE_JSON) as json_file:
        attribute_dictionary = json.load(json_file)
    variable_attribute_dictionary = attribute_dictionary['_variables']
    global_attribute_dictionary = attribute_dictionary['_global']

    ## set variable attrs
    for var in list(ds.variables):
        ds[var].setncatts(variable_attribute_dictionary[var])

    ## set global attrs
    timeformat = '%Y-%m-%dT%H:%M:%SZ'
    file_timeformat = '%Y%m%d'

    time_start = nc4.num2date(np.min(TIME[:]), time_units, time_calendar).strftime(timeformat)
    time_end = nc4.num2date(np.max(TIME[:]), time_units, time_calendar).strftime(timeformat)
    time_start_filename = nc4.num2date(np.min(TIME[:]), time_units, time_calendar).strftime(file_timeformat)
    time_end_filename = nc4.num2date(np.max(TIME[:]), time_units, time_calendar).strftime(file_timeformat)

    contributor_name, contributor_email, contributor_role = get_contributors(files_to_agg)
    add_attribute = {
                    'title':                    ("Long Timeseries Velocity Aggregated product: " + VoI + " at " +
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
                    'keywords':                 ', '.join([VoI, 'AGGREGATED']),
                    'rejected_files':           "\n".join(rejected_files),
                    'contributor_name':         "; ".join(contributor_name),
                    'contributor_email':        "; ".join(contributor_email),
                    'contributor_role':         "; ".join(contributor_role)}

    github_comment = ('\nThis file was created using https://github.com/aodn/python-aodntools/blob/'
                      '{v}/aodntools/timeseries_products/aggregated_timeseries.py'.format(v=__version__)
                      )
    global_attribute_dictionary['lineage'] += github_comment

    global_attribute_dictionary.update(add_attribute)
    ds.setncatts(dict(sorted(global_attribute_dictionary.items())))

    ds.close()


    ## create the output file name and rename the tmp file
    facility_code = get_facility_code(files_to_agg[0])
    data_code = get_data_code(VoI)
    product_type = 'aggregated-timeseries'
    file_version = 1
    output_name = '_'.join(['IMOS', facility_code, data_code, time_start_filename, site_code, ('FV0'+str(file_version)),
                            (VoI+"-"+product_type),
                            ('END-'+ time_end_filename), 'C-' + datetime.utcnow().strftime(file_timeformat)]) + '.nc'
    ncout_path = os.path.join(base_path, output_name)
    os.rename(outfile, ncout_path)

    return ncout_path, bad_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenate X,Y,Z velocity variables from ALL instruments from ALL deployments from ONE site")
    parser.add_argument('-site', dest='site_code', help='site code, like NRMMAI',  required=True)
    parser.add_argument('-var', dest='VoI', help='variable to aggregate, like TEMP', required=True)
    parser.add_argument('-files', dest='filenames', help='name of the file that contains the source URLs', required=True)
    parser.add_argument('-path', dest='output_path', help='path where the result file will be written. Default: ./', default='./', required=False)
    args = parser.parse_args()

    with open(args.filenames) as ff:
        files_to_agg = [line.rstrip() for line in ff]


    print(aggregate_timeseries(files_to_agg=files_to_agg, site_code=args.site_code, VoI=args.VoI, base_path = args.output_path))
