import os
import tempfile
import shutil
from copy import deepcopy

from  netCDF4 import Dataset, num2date, stringtochar
import numpy as np
import json
import argparse
from pkg_resources import resource_filename
from aodntools import __version__

import xarray as xr

from aodntools.timeseries_products import aggregated_timeseries as utils
from aodntools.timeseries_products.common import (NoInputFilesError, check_velocity_file, current_utc_timestamp,
                                                  TIMESTAMP_FORMAT, DATESTAMP_FORMAT)

TEMPLATE_JSON = resource_filename(__name__, 'velocity_aggregated_timeseries_template.json')


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
    :return: file path of the aggregated product, dict of rejected files: errors
    """

    varlist = ['UCUR', 'VCUR', 'WCUR', 'DEPTH']
    time_units="days since 1950-01-01 00:00:00 UTC"
    time_calendar="gregorian"
    epoch = np.datetime64("1950-01-01T00:00:00")
    one_day = np.timedelta64(1, 'D')

    bad_files = {}

    # default name for temporary file. It will be renamed at the end
    fd, temp_outfile = tempfile.mkstemp(suffix='.nc', dir=output_dir)

    ## check files and get total number of flattened obs
    n_obs_total = 0
    for file in files_to_agg:
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            error_list = check_velocity_file(nc, site_code)
            if not error_list:
                nc = utils.in_water(nc)
                n_obs_total += get_number_flatvalues(nc)[0]
            else:
                bad_files.update({file: error_list})

    # remove bad files form the list and sort in chronological order
    files_to_agg = deepcopy(files_to_agg)
    for file in bad_files.keys():
        files_to_agg.remove(file)
    if len(files_to_agg) == 0:
        raise NoInputFilesError("no valid input files to aggregate")
    files_to_agg = utils.sort_files(files_to_agg, input_dir=input_dir)

    n_files = len(files_to_agg)

    ## create ncdf file, dimensions and variables
    ds = Dataset(os.path.join(output_dir, temp_outfile), 'w', format="NETCDF4_CLASSIC")
    OBSERVATION = ds.createDimension('OBSERVATION', size=n_obs_total)
    INSTRUMENT = ds.createDimension('INSTRUMENT', size=n_files)
    STRING256 = ds.createDimension("strlen", 256)

    obs_double_template = {'datatype': np.float64, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_float_template = {'datatype': np.float32, 'zlib': True, 'dimensions': ('OBSERVATION'), "fill_value": 99999.0}
    obs_byte_template = {'datatype': np.byte, 'zlib': True, 'dimensions': ('OBSERVATION'), 'fill_value': 99}
    obs_int_template = {'datatype': np.int16, 'zlib': True, 'dimensions': ('OBSERVATION')}
    inst_S256_template = {'datatype': 'S1', 'dimensions': ('INSTRUMENT', "strlen")}
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
    CELL_INDEX = ds.createVariable(varname='CELL_INDEX', **obs_int_template)

    ## main loop
    start = 0
    for index, file in enumerate(files_to_agg):
        print(index)
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:

            ## in water data only
            nc = utils.in_water(nc)
            n_measurements = len(nc.TIME)

            n_obs, n_cells = get_number_flatvalues(nc)
            end = start + n_obs
            UCUR[start:end] = flat_variable(nc, 'UCUR')
            UCURqc[start:end] = flat_variable(nc, 'UCUR_quality_control')
            VCUR[start:end] = flat_variable(nc, 'VCUR')
            VCURqc[start:end] = flat_variable(nc, 'VCUR_quality_control')
            if 'WCUR' in nc.data_vars:
                WCUR[start:end] = flat_variable(nc, 'WCUR')
                WCURqc[start:end] = flat_variable(nc, 'WCUR_quality_control')
            else:
                WCUR[start:end] = np.ma.masked
                WCURqc[start:end] = np.full(n_obs, 9)

            ##calculate depth and add CELL_INDEX
            if 'HEIGHT_ABOVE_SENSOR' in nc.dims:
                DEPTH[start:end] = (nc.DEPTH - nc.HEIGHT_ABOVE_SENSOR).values.flatten()
                DEPTHqc[start:end] = np.repeat(flat_variable(nc, 'DEPTH_quality_control'), n_cells)
                CELL_INDEX[start:end] = np.tile(np.arange(n_cells, dtype=np.uint32), n_measurements)

            else:
                DEPTH[start:end] = nc.DEPTH.values
                DEPTHqc[start:end] = nc.DEPTH_quality_control.values
                CELL_INDEX[start:end] = np.full(n_obs, 0, dtype=np.uint32)

            ## set TIME and instrument index
            TIME[start:end] = (np.repeat(flat_variable(nc, 'TIME'), n_cells) - epoch) / one_day
            instrument_index[start:end] = np.repeat(index, n_obs)
            ## get and store deployment metadata
            LATITUDE[index] = nc.LATITUDE.values
            LONGITUDE[index] = nc.LONGITUDE.values
            NOMINAL_DEPTH[index] = utils.get_nominal_depth(nc)
            source_file[index] = stringtochar(np.array(file, dtype='S256'))
            instrument_id[index] = stringtochar(np.array(utils.get_instrument_id(nc), dtype='S256'))
            ## add time offset to the middle of the measuring window, if it exists
            if 'seconds_to_middle_of_measurement' in nc.TIME.attrs:
                SECONDS_TO_MIDDLE[index] = nc.TIME.seconds_to_middle_of_measurement
            else:
                SECONDS_TO_MIDDLE[index] = np.nan
        start = end

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
    time_start = num2date(np.min(TIME[:]), time_units, time_calendar).strftime(TIMESTAMP_FORMAT)
    time_end = num2date(np.max(TIME[:]), time_units, time_calendar).strftime(TIMESTAMP_FORMAT)
    time_start_filename = num2date(np.min(TIME[:]), time_units, time_calendar).strftime(DATESTAMP_FORMAT)
    time_end_filename = num2date(np.max(TIME[:]), time_units, time_calendar).strftime(DATESTAMP_FORMAT)

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
                    'date_created':             current_utc_timestamp(),
                    'history':                  current_utc_timestamp() + ': Aggregated file created.',
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

    ds.close()


    ## create the output file name and rename the tmp file
    facility_code = utils.get_facility_code(os.path.join(input_dir, files_to_agg[0]))
    data_code = 'VZ'
    product_type = 'aggregated-timeseries'
    file_version = 1
    output_name = '_'.join(['IMOS', facility_code, data_code, time_start_filename, site_code, ('FV0'+str(file_version)),
                            ("velocity-"+product_type),
                            ('END-'+ time_end_filename), 'C-' + current_utc_timestamp(DATESTAMP_FORMAT)]) + '.nc'
    ncout_path = os.path.join(output_dir, output_name)
    os.close(fd)
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

    print(velocity_aggregated(files_to_agg=files_to_agg, site_code=args.site_code,
                              input_dir=args.input_dir, output_dir=args.output_dir,
                              download_url_prefix=args.download_url, opendap_url_prefix=args.opendap_url))
