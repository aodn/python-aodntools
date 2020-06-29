import argparse
import json
import os
import shutil
import tempfile
from datetime import datetime

import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset, num2date, stringtochar
from pkg_resources import resource_filename

import aodntools.timeseries_products.aggregated_timeseries as utils
from aodntools import __version__
from aodntools.timeseries_products.common import NoInputFilesError
from aodntools.timeseries_products.velocity_aggregated_timeseries import check_file

TEMPLATE_JSON = resource_filename(__name__, 'velocity_hourly_timeseries_template.json')
QC_FLAG_MAX = 2
TIME_UNITS = "days since 1950-01-01 00:00:00 UTC"
TIME_CALENDAR = "gregorian"
TIME_EPOCH = np.datetime64("1950-01-01T00:00:00")
ONE_DAY = np.timedelta64(1, 'D')


def cell_velocity_resample(df, binning_function):
    """
    Resample a dataset to a specific time_interval.
    if WCUR not present, returns nan
    :param df: grouped dataframe
    :param binning_function: name of standard numpy function used for binning
    :return: binned U, v, W CUR according to the binning function
    """
    df_binned = df.apply(binning_function)
    UCUR = np.array(df_binned['UCUR'])
    VCUR = np.array(df_binned['VCUR'])
    if 'WCUR' in df_binned:
        WCUR = np.array(df_binned['WCUR'])
    else:
        WCUR = np.full(len(df), np.nan)
    DEPTH = np.array(df_binned['DEPTH'])

    return UCUR, VCUR, WCUR, DEPTH


def append_resampled_values(nc_cell, ds, slice_start, binning_functions):
    """
    Resample U, V, W current and depth values from a single ADCP cell into hourly bins, and
    append the mean values to the corresponding variables in the output dataset (starting at
    index slice_start), along with additional statistical variables specified by binning_functions.
    :param nc_cell: input xarray Dataset representing a single ADCP cell (or point time series)
    :param ds: output netcdf4 Dataset to update with resampled values
    :param slice_start: start index of the slice
    :param binning_functions: list of numpy function names for binning
    :return: end index of the slice
    """
    df_cell = nc_cell.squeeze().to_dataframe()
    # shift the index forward 30min to centre the bins on the hour
    df_cell.index = df_cell.index + pd.Timedelta(minutes=30)
    # TODO: shift timestamps to centre of sampling interval

    df_cell_1H = df_cell.resample('1H')
    slice_end = len(df_cell_1H) + slice_start

    # set binned timestamps
    time_slice = (np.fromiter(df_cell_1H.groups.keys(), dtype='M8[ns]') - TIME_EPOCH) / ONE_DAY
    ds['TIME'][slice_start:slice_end] = time_slice

    # take the mean of the variables
    ds['UCUR'][slice_start:slice_end], \
    ds['VCUR'][slice_start:slice_end], \
    ds['WCUR'][slice_start:slice_end], \
    ds['DEPTH'][slice_start:slice_end] = cell_velocity_resample(df_cell_1H, 'mean')

    for method in binning_functions:
        ds['UCUR_' + method][slice_start:slice_end], \
        ds['VCUR_' + method][slice_start:slice_end], \
        ds['WCUR_' + method][slice_start:slice_end], \
        ds['DEPTH_' + method][slice_start:slice_end] = cell_velocity_resample(df_cell_1H, method)

    return slice_end



## MAIN FUNCTION
def velocity_hourly_aggregated(files_to_agg, site_code, input_dir='', output_dir='./',
                               download_url_prefix=None, opendap_url_prefix=None):
    """
    Aggregate U, V and W CUR variables from the given files (from the same site) and average into hourly bins.
    The vertical cells are flattened and the actual depth of each is calculated.
    Additional metadata variables are stored to track the origin of the data.
    :param files_to_agg: list of files to aggregate
    :param site_code: site code
    :param input_dir: base path where source files are stored
    :param output_dir: path where the result file will be written
    :param download_url_prefix: URL prefix for file download (to be prepended to paths in files_to_agg)
    :param opendap_url_prefix: URL prefix for OPENAP access (to be prepended to paths in files_to_agg)
    :return: file path of the hourly aggregated product, dict of rejected files: errors
    """

    varlist = ['UCUR', 'VCUR', 'WCUR', 'DEPTH']
    binning_fun = ['max', 'min', 'std', 'count']

    bad_files = {}

    chunk_size = 90  ## size in days

    ## default name for temporary file. It will be renamed at the end
    _, temp_outfile = tempfile.mkstemp(suffix='.nc', dir=output_dir)

    ## check files and get total number of flattened obs
    print("CHECKING FILES...")
    for index, file in enumerate(files_to_agg):
        print(index, end=',', flush=True)
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            error_list = check_file(nc, site_code)
            if error_list:
                bad_files.update({file: error_list})
    print(" ")

    ## remove bad files form the list
    for file in bad_files.keys():
        files_to_agg.remove(file)
    if len(files_to_agg) == 0:
        raise NoInputFilesError("no valid input files to aggregate")

    ## sort the files in chronological order
    files_to_agg = utils.sort_files(files_to_agg, input_dir=input_dir)

    ## create ncdf file, dimensions (unlimited) and variables
    ds = Dataset(temp_outfile, 'w', format='NETCDF4_CLASSIC')
    OBSERVATION = ds.createDimension('OBSERVATION', size=None)
    INSTRUMENT = ds.createDimension('INSTRUMENT', size=len(files_to_agg))
    STRING256 = ds.createDimension("strlen", size=256)

    obs_double_template = {'datatype': np.float64, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_float_template = {'datatype': np.float32, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_int_template = {'datatype': np.int16, 'zlib': True, 'dimensions': ('OBSERVATION')}
    inst_S256_template = {'datatype': 'S1', 'dimensions': ('INSTRUMENT', "strlen")}
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
    SECONDS_TO_MIDDLE = ds.createVariable(varname='SECONDS_TO_MIDDLE', **inst_float_template)
    CELL_INDEX = ds.createVariable(varname='CELL_INDEX', **obs_int_template)



    ## main loop
    print('PROCESSING...')
    slice_start = 0
    for index, file in enumerate(files_to_agg):
        print(index, end=",", flush=True)

        ## this is for filling the slice of variables with INSTRUMENT dim
        slice_instrument_start = slice_start

        with xr.open_dataset(os.path.join(input_dir, file)) as nc:

            is_2D = 'HEIGHT_ABOVE_SENSOR' in list(nc.variables)

            ## mask values with QC flag>2
            for var in varlist:
                nc[var] = nc[var].where(nc[var+'_quality_control'] <= QC_FLAG_MAX)

            ## process in chunks
            ## in water only
            chunk_start = np.datetime64(nc.attrs['time_deployment_start'])
            chunk_end = np.datetime64(nc.attrs['time_deployment_end'])

            time_increment = 60*60*24*chunk_size    ## secs x mins x hours x days
            chunk_increment = np.timedelta64(time_increment, 's')
            chunk_partial = chunk_start + chunk_increment
            chunk_index = 0
            while chunk_start < chunk_partial and chunk_start <= chunk_end:
                nc_chunk = nc.where((nc.TIME >= chunk_start) & (nc.TIME < chunk_partial), drop=True)
                if is_2D:
                    ## process all cells, one by one
                    heights = nc_chunk.HEIGHT_ABOVE_SENSOR.values
                    for cell_idx, cell_height in enumerate(heights):
                        ## get cell data, drop HEIGHT_ABOVE_SENSOR dim
                        nc_cell = nc_chunk.sel(HEIGHT_ABOVE_SENSOR=cell_height)
                        ## convert to absolute DEPTH
                        nc_cell['DEPTH'] = nc_cell['DEPTH'] - cell_height
                        slice_end = append_resampled_values(nc_cell[varlist], ds, slice_start, binning_fun)
                        CELL_INDEX[slice_start:slice_end] = np.full(slice_end - slice_start, cell_idx, dtype=np.uint32)

                        slice_start = slice_end
                else:
                    slice_end = append_resampled_values(nc_chunk[varlist], ds, slice_start, binning_fun)
                    CELL_INDEX[slice_start:slice_end] = np.full(slice_end - slice_start, 0, dtype=np.uint32)

                    slice_start = slice_end
                chunk_start = chunk_partial
                chunk_partial += chunk_increment
                chunk_index += 1


            ## metadata variables
            instrument_index[slice_instrument_start:slice_end] = np.repeat(index, slice_end - slice_instrument_start)
            LATITUDE[index] = nc.LATITUDE.values
            LONGITUDE[index] = nc.LONGITUDE.values
            NOMINAL_DEPTH[index] = np.array(utils.get_nominal_depth(nc))
            source_file[index] = stringtochar(np.array(file, dtype='S256'))
            instrument_id[index] = stringtochar(np.array(utils.get_instrument_id(nc), dtype='S256'))
            ## add time offset to the middle of the measuring window, if it exists
            if 'seconds_to_middle_of_measurement' in nc.TIME.attrs:
                SECONDS_TO_MIDDLE[index] = nc.TIME.seconds_to_middle_of_measurement
            else:
                SECONDS_TO_MIDDLE[index] = np.nan

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

    time_start = num2date(np.min(TIME[:]), TIME_UNITS, TIME_CALENDAR).strftime(timeformat)
    time_end = num2date(np.max(TIME[:]), TIME_UNITS, TIME_CALENDAR).strftime(timeformat)
    time_start_filename = num2date(np.min(TIME[:]), TIME_UNITS, TIME_CALENDAR).strftime(file_timeformat)
    time_end_filename = num2date(np.max(TIME[:]), TIME_UNITS, TIME_CALENDAR).strftime(file_timeformat)


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
                    'rejected_files':           "\n".join(bad_files.keys()),
                    'generating_code_version':  __version__
    }
    add_attribute.update(utils.get_contributors(files_to_agg=files_to_agg, input_dir=input_dir))

    ## add version
    github_comment = ('\nThis file was created using https://github.com/aodn/python-aodntools/blob/'
                      '{v}/aodntools/timeseries_products/{f}'.format(v=__version__, f=os.path.basename(__file__))
                      )
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
    file_version = 2
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


    print(velocity_hourly_aggregated(files_to_agg=files_to_agg, site_code=args.site_code,
                                     input_dir=args.input_dir, output_dir=args.output_dir,
                                     download_url_prefix=args.download_url, opendap_url_prefix=args.opendap_url))
