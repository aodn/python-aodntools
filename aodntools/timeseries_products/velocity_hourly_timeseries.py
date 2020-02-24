import os
import sys
import tempfile
import shutil
from netCDF4 import Dataset, num2date
import numpy as np
import json
from datetime import datetime
import argparse

from pkg_resources import resource_filename
from aodntools import __version__
import aodntools.timeseries_products.aggregated_timeseries as utils

import xarray as xr
import pandas as pd


TEMPLATE_JSON = resource_filename(__name__, 'velocity_hourly_timeseries_template.json')


def check_file(nc, site_code):
    """
    Return list of errors found in the file:
    Variables of interest are present
    TIME, DEPTH, LATITUDE, LONGITUDE,  is present
    NOMINAL_DEPTH is not present as variable or attribute
    file_version is not FV01
    if LATITUDE or LONIGITUDE dimension has length >1
    if TIME.seconds_to_middle_of_measurement exist

    :param nc: xarray dataset
    :param site_code: code of the mooring site
    :return: dictionary with the file name and list of failed tests
    """

    attributes = list(nc.attrs)
    variables = list(nc.variables)
    allowed_dimensions = ['TIME', 'LATITUDE', 'LONGITUDE', 'HEIGHT_ABOVE_SENSOR']
    required_variables = ['UCUR', 'VCUR', 'WCUR', 'DEPTH']
    error_list = []

    if nc.site_code != site_code:
        error_list.append('Wrong site_code: ' + nc.site_code)

    if 'Level 1' not in nc.file_version:
        error_list.append('Wrong file version: ' + nc.file_version)
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
    if 'seconds_to_middle_of_measurement' not in nc['TIME'].attrs:
        error_list.append('seconds_to_middle_of_measurement not present in TIME')

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

def cell_velocity_resample(df, binning_function, is_WCUR):
    """
    Resample a dataset to a specific time_interval.
    if WCUR not present, returns nan
    :param df: grouped dataframe
    :param binning_function: function used for binning as non string standard numpy function
    :param is_WCUR: True if WCUR is present in nc, False otherwise
    :return: binned U, v, W CUR according to the binning function
    """
    df_binned = df.apply(binning_function)
    UCUR = df_binned['UCUR']
    VCUR = df_binned['VCUR']
    if is_WCUR:
        WCUR = df_binned['WCUR']
    else:
        WCUR = np.full(len(df), np.nan)
    DEPTH = df_binned['DEPTH']

    return UCUR, VCUR, WCUR, DEPTH


def get_resampled_values(nc_cell, ds, slice_start, varlist, binning_fun, epoch, one_day, is_WCUR):
    """
    get U, V, W current values resampled
    :param nc_cell: xarray DATASET
    :param ds: netcdf4 dataset
    :param slice_start: start index of the slice
    :param varlist: list of variable names to subset the dataset
    :param binnig_fun: list of function names for binning
    :param one_day: timedelta one day
    :param epoch: base epoch
    :param is_WCUR: flag indicating if WCUR is present
    :return: end index of the slice
    """
    nc_cell = nc_cell.where(nc_cell.DEPTH_quality_control < 4, drop=True)
    nc_cell = nc_cell[varlist]
    nc_cell = nc_cell.to_dataframe()
    ## back the index 30min
    nc_cell.index = nc_cell.index - pd.Timedelta(30, units='m')

    nc_cell_1H = nc_cell.resample('1H')
    slice_end = len(nc_cell_1H) + slice_start

    ## move time it forward and get it
    time_slice = ((np.fromiter(nc_cell_1H.groups.keys(), dtype='M8[ns]') + np.timedelta64(1, 'h')) - epoch) / one_day
    ds['TIME'][slice_start:slice_end] = time_slice

    # take the mean of the variables
    ds['UCUR'][slice_start:slice_end], \
    ds['VCUR'][slice_start:slice_end], \
    ds['WCUR'][slice_start:slice_end], \
    ds['DEPTH'][slice_start:slice_end] = cell_velocity_resample(nc_cell_1H, 'mean', is_WCUR)

    for method in binning_fun:
        ds['UCUR_' + method][slice_start:slice_end], \
        ds['VCUR_' + method][slice_start:slice_end], \
        ds['WCUR_' + method][slice_start:slice_end], \
        ds['DEPTH_' + method][slice_start:slice_end] = cell_velocity_resample(nc_cell_1H, method, is_WCUR)

    return slice_end



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
    binning_fun = ['max', 'min', 'std', 'count']
    
    time_units="days since 1950-01-01 00:00:00 UTC"
    time_calendar="gregorian"
    epoch = np.datetime64("1950-01-01T00:00:00")
    one_day = np.timedelta64(1, 'D')

    bad_files = []
    rejected_files = []

    ## default name for temporary file. It will be renamed at the end
    _, temp_outfile = tempfile.mkstemp(suffix='.nc', dir=output_dir)

    ## check files and get total number of flattened obs
    print("CHECKING FILES...")
    for index, file in enumerate(files_to_agg):
        print(index, end=',', flush=True)
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            nc = utils.in_water(nc)
            error_list = check_file(nc, site_code)
            if error_list:
                bad_files.append([file, error_list])
                rejected_files.append(file)
    print(" ")

    ## remove bad files form the list
    for file in bad_files:
        files_to_agg.remove(file[0])

    ## sort the files in chronological order
    files_to_agg = utils.sort_files(files_to_agg, input_dir=input_dir)


    ## create ncdf file, dimensions (unlimited) and variables
    ds = Dataset(os.path.join(output_dir, temp_outfile), 'w')
    OBSERVATION = ds.createDimension('OBSERVATION', size=None)
    INSTRUMENT = ds.createDimension('INSTRUMENT', size=None)

    obs_double_template = {'datatype': np.float64, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_float_template = {'datatype': np.float32, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_int_template = {'datatype': np.uint32, 'zlib': True, 'dimensions': ('OBSERVATION')}
    inst_S256_template = {'datatype': 'str', 'dimensions': ('INSTRUMENT')}
    inst_float_template ={'datatype': np.float32, 'dimensions': ('INSTRUMENT')}
    inst_double_template ={'datatype': np.float64, 'dimensions': ('INSTRUMENT')}

    UCUR = ds.createVariable(varname='UCUR', **obs_float_template)
    UCUR_max = ds.createVariable(varname='UCUR_max', **obs_float_template)
    UCUR_min = ds.createVariable(varname='UCUR_min', **obs_float_template)
    UCUR_std = ds.createVariable(varname='UCUR_std', **obs_float_template)
    UCUR_count = ds.createVariable(varname='UCUR_count', **obs_int_template)
    VCUR = ds.createVariable(varname='VCUR', **obs_float_template)
    VCUR_max = ds.createVariable(varname='VCUR_max', **obs_float_template)
    VCUR_min = ds.createVariable(varname='VCUR_min', **obs_float_template)
    VCUR_std = ds.createVariable(varname='VCUR_std', **obs_float_template)
    VCUR_count = ds.createVariable(varname='VCUR_count', **obs_int_template)
    WCUR = ds.createVariable(varname='WCUR', **obs_float_template)
    WCUR_max = ds.createVariable(varname='WCUR_max', **obs_float_template)
    WCUR_min = ds.createVariable(varname='WCUR_min', **obs_float_template)
    WCUR_std = ds.createVariable(varname='WCUR_std', **obs_float_template)
    WCUR_count = ds.createVariable(varname='WCUR_count', **obs_int_template)

    DEPTH = ds.createVariable(varname='DEPTH', **obs_float_template)
    DEPTH_max = ds.createVariable(varname='DEPTH_max', **obs_float_template)
    DEPTH_min = ds.createVariable(varname='DEPTH_min', **obs_float_template)
    DEPTH_std = ds.createVariable(varname='DEPTH_std', **obs_float_template)
    DEPTH_count = ds.createVariable(varname='DEPTH_count', **obs_int_template)

    TIME = ds.createVariable(varname='TIME', **obs_double_template)
    instrument_index = ds.createVariable(varname='instrument_index', **obs_int_template)

    source_file = ds.createVariable(varname='source_file', **inst_S256_template)
    instrument_id = ds.createVariable(varname='instrument_id', **inst_S256_template)
    LATITUDE = ds.createVariable(varname='LATITUDE', **inst_double_template)
    LONGITUDE = ds.createVariable(varname='LONGITUDE', **inst_double_template)
    NOMINAL_DEPTH = ds.createVariable(varname='NOMINAL_DEPTH', **inst_float_template)
    CELL_INDEX = ds.createVariable(varname='CELL_INDEX', **obs_int_template)

    #SECONDS_TO_MIDDLE = ds.createVariable(varname='SECONDS_TO_MIDDLE', **inst_float_template)


    ## main loop
    print('PROCESSING...')
    slice_start = 0
    for index, file in enumerate(files_to_agg):
        print(index, end=",", flush=True)

        ## this is for filling the slice of variables with INSTRUMENT dim
        slice_instrument_start = slice_start

        with xr.open_dataset(os.path.join(input_dir, file)) as nc:

            if 'HEIGHT_ABOVE_SENSOR' in list(nc.variables):
                is_2D = True
            else:
                is_2D = False

            if 'WCUR' in list(nc.data_vars):
                is_WCUR = True
            else:
                is_WCUR = False

            ## in-water data only
            nc = utils.in_water(nc)
            n_measurements = len(nc.TIME)

            ## move timestamp to the middle of the measurement window
            time_delta_ns = int(nc['TIME'].seconds_to_middle_of_measurement * 10**9)
            nc['TIME'] = nc['TIME'] + np.timedelta64(time_delta_ns, 'ns')

            if is_2D:
                ## process all cells, one by one
                cells = nc.HEIGHT_ABOVE_SENSOR.values
                for cell_idx, cell in enumerate(cells):
                    ## get cell data, drop HEIGHT_ABOVE_SENSOR dim
                    nc_cell = nc.where(nc.HEIGHT_ABOVE_SENSOR == cell, drop=True).squeeze('HEIGHT_ABOVE_SENSOR')
                    ## convert to absolute DEPTH
                    nc_cell['DEPTH'] = nc_cell['DEPTH'] - cell
                    slice_end = get_resampled_values(nc_cell, ds, slice_start, varlist, binning_fun,
                                                     epoch, one_day, is_WCUR)
                    CELL_INDEX[slice_start:slice_end] = np.full(slice_end - slice_start, cell_idx, dtype=np.uint32)

                    slice_start = slice_end
            else:
                slice_end = get_resampled_values(nc_cell, ds, slice_start, varlist, binning_fun,
                                                 epoch, one_day, is_WCUR)
                CELL_INDEX[slice_start:slice_end] = np.full(slice_end - slice_start, 0, dtype=np.uint32)

                slice_start = slice_end

            ## metadata variables
            instrument_index[slice_instrument_start:slice_end] = np.repeat(index, slice_end - slice_instrument_start)
            LATITUDE[index] = nc.LATITUDE.values
            LONGITUDE[index] = nc.LONGITUDE.values
            NOMINAL_DEPTH[index] = np.array(utils.get_nominal_depth(nc))
            instrument_id[index] = utils.get_instrument_id(nc)
            source_file[index] = file

    print(" ")
    ## add atributes
    with open(TEMPLATE_JSON) as json_file:
        attribute_dictionary = json.load(json_file)
    variable_attribute_dictionary = attribute_dictionary['_variables']
    global_attribute_dictionary = attribute_dictionary['_global']

    ## set variable attrs
    for var in list(ds.variables):
        ds[var].setncatts(variable_attribute_dictionary[var])

    if download_url_prefix or opendap_url_prefix:
        ds['source_file'].setncatts(utils.source_file_attributes(download_url_prefix, opendap_url_prefix))

    ## set global attrs
    timeformat = '%Y-%m-%dT%H:%M:%SZ'
    file_timeformat = '%Y%m%d'

    time_start = num2date(np.min(TIME[:]), time_units, time_calendar).strftime(timeformat)
    time_end = num2date(np.max(TIME[:]), time_units, time_calendar).strftime(timeformat)
    time_start_filename = num2date(np.min(TIME[:]), time_units, time_calendar).strftime(file_timeformat)
    time_end_filename = num2date(np.max(TIME[:]), time_units, time_calendar).strftime(file_timeformat)


    contributor_name, contributor_email, contributor_role = utils.get_contributors(files_to_agg=files_to_agg, input_dir=input_dir)
    add_attribute = {
                    'title':                    ("Long Timeseries Velocity Hourly Aggregated product: " + ', '.join(varlist) + " at " +
                                                  site_code + " between " + time_start + " and " + time_end),
                    'site_code':                site_code,
                    'time_coverage_start':      time_start,
                    'time_coverage_end':        time_end,
                    'geospatial_vertical_min':  np.float32(np.nanmin(ds['DEPTH'])),
                    'geospatial_vertical_max':  np.float32(np.nanmax(ds['DEPTH'])),
                    'geospatial_lat_min':       np.float64(np.min(ds['LATITUDE'])),
                    'geospatial_lat_max':       np.float64(np.max(ds['LATITUDE'])),
                    'geospatial_lon_min':       np.float64(np.min(ds['LONGITUDE'])),
                    'geospatial_lon_max':       np.float64(np.max(ds['LONGITUDE'])),
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
                      '{v}/aodntools/timeseries_products/aggregated_timeseries.py'.format(v=__version__))
    global_attribute_dictionary['lineage'] += github_comment

    global_attribute_dictionary.update(add_attribute)
    ds.setncatts(dict(sorted(global_attribute_dictionary.items())))


    ## NOTE: There is a possibility of having NaNs in DEPTH after the binning
    ## this is the warning when calculating the min/max DEPTH
    ## maybe I should clean the dataset before close it

    ds.close()



    ## create the output file name and rename the tmp file
    facility_code = utils.get_facility_code(os.path.join(input_dir, files_to_agg[0]))
    data_code = 'VZ'
    product_type = 'hourly-timeseries'
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
