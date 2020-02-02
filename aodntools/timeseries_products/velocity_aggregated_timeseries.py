import os
import sys
import netCDF4 as nc4
import numpy as np
import json
from datetime import datetime
import argparse

import xarray as xr

import aggregated_timeseries as TStools

TEMPLATE_JSON = 'velocity_aggregated_timeseries_template.json'


def sort_files(files_to_agg):
    """
    sort list of files according to deployment date
    requires netcdf4 as nc4, dateutil.parser as parse
    :param files_to_agg: List of files to sort
    :return: sorted list of files
    """

    time_start = []
    for file in files_to_agg:
        with nc4.Dataset(file, 'r') as ds:
            time_start.append(np.datetime64(ds.time_deployment_start))
    tuples = sorted(zip(time_start, files_to_agg))
    return [t[1] for t in tuples]


def check_file(nc, site_code):
    """
    Return list of errors found in the file:
    Variables of interest are present
    TIME. LATITUDE, LONGITUDE,  is present
    NOMINAL_DEPTH is not present as variable or attribute
    file_version is not FV01
    if LATITUDE or LONIGUTDE dimension has length >1

    :param file: name of the netcdf file
    :param nc: xarray dataset
    :param VoI: string. Variable of Interest
    :param site_code: code of the mooring site
    :return: dictionary with the file name and list of failed tests
    """

    attributes = list(nc.attrs)
    file_variables = list(nc.variables)
    allowed_dimensions = ['TIME', 'LATITUDE', 'LONGITUDE', 'HEIGHT_ABOVE_SENSOR']
    required_variables = ['UCUR', 'VCUR', 'WCUR']
    error_list = []

    nc_site_code = nc.site_code
    if nc_site_code != site_code:
        error_list.append('Wrong site_code: ' + nc_site_code)

    nc_file_version = nc.file_version
    if 'Level 1' not in nc_file_version:
        error_list.append('Wrong file version: ' + nc_file_version)

    for variable in required_variables:
        if variable not in file_variables:
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

    if 'NOMINAL_DEPTH' not in file_variables and 'instrument_nominal_depth' not in attributes:
        error_list.append('no NOMINAL_DEPTH')

    return error_list


def get_nvalues(nc):
    """
    Get the number of cells above the sensor
    :param nc: xarray dataset
    :return: number of cells, number of bins
    """
    if 'HEIGHT_ABOVE_SENSOR' in nc.dims:
        nbins = nc.dims['HEIGHT_ABOVE_SENSOR']
        nvalues = nc.dims['TIME'] * nbins
    else:
        nbins = 1
        nvalues = nc.dims['TIME']
    return nvalues, nbins


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
    cut data to in-water only timestamps, dropping the out-of-water records.
    :param nc: xarray dataset
    :return: xarray dataset
    """
    time_deployment_start = np.datetime64(nc.attrs['time_deployment_start'][:-1])
    time_deployment_end = np.datetime64(nc.attrs['time_deployment_end'][:-1])
    TIME = nc['TIME'][:]
    return nc.where((TIME >= time_deployment_start) & (TIME <= time_deployment_end), drop=True)




## MAIN FUNCTION
def aggregate_velocity(files_to_agg, site_code, base_path):
    """
    Aggregate U, V and W CUR variables from all deployments at one site.
    the vertical cells are flattened and related to its depth
    additional metadata variables are stored to track the origin of the data
    :param files_to_agg: list of files to aggregate
    :param site_code: site code
    :param base_path: base path to store the resulting file
    :return: name of the resulting file, list of rejected files
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
    outfile = 'Velocity_agg_tmp.nc'

    ## sort the file list in chronological order
    files_to_agg = sort_files(files_to_agg)

    ## check files and get total number of flattened obs
    for file in files_to_agg:
        with xr.open_dataset(file) as nc:
            ## clip to in water data only
            nc = in_water(nc)

            varlen_file.append(get_nvalues(nc))
            error_list = check_file(nc, site_code)
            if not error_list:
                varlen_list.append(get_nvalues(nc)[0])
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
    ds = nc4.Dataset(os.path.join(base_path, outfile), 'w')
    OBSERVATION = ds.createDimension('OBSERVATION', size=varlen_total)
    INSTRUMENT = ds.createDimension('INSTRUMENT', size=n_files)

    UCUR = ds.createVariable(varname='UCUR', datatype='float', zlib=True, dimensions=('OBSERVATION'), fill_value=999999.0)
    VCUR = ds.createVariable(varname='VCUR', datatype='float', zlib=True, dimensions=('OBSERVATION'), fill_value=999999.0)
    WCUR = ds.createVariable(varname='WCUR', datatype='float', zlib=True, dimensions=('OBSERVATION'), fill_value=999999.0)
    DEPTH = ds.createVariable(varname='DEPTH', datatype='float', zlib=True, dimensions=('OBSERVATION'), fill_value=999999.0)
    UCURqc = ds.createVariable(varname='UCUR_quality_control', datatype='byte', zlib=True, dimensions=('OBSERVATION'), fill_value=99)
    VCURqc = ds.createVariable(varname='VCUR_quality_control', datatype='byte', zlib=True, dimensions=('OBSERVATION'), fill_value=99)
    WCURqc = ds.createVariable(varname='WCUR_quality_control', datatype='byte', zlib=True, dimensions=('OBSERVATION'), fill_value=99)
    DEPTHqc = ds.createVariable(varname='DEPTH_quality_control', datatype='byte', zlib=True, dimensions=('OBSERVATION'), fill_value=99)
    TIME = ds.createVariable(varname='TIME', datatype='double', zlib=True, dimensions=('OBSERVATION'))
    instrument_index = ds.createVariable(varname='instrument_index', datatype='int', zlib=True, dimensions=('OBSERVATION'))

    source_file = ds.createVariable(varname='source_file', datatype='S256', dimensions=('INSTRUMENT'))
    instrument_id = ds.createVariable(varname='instrument_id', datatype='S256', dimensions=('INSTRUMENT'))
    LATITUDE = ds.createVariable(varname='LATITUDE', datatype='float', dimensions=("INSTRUMENT"))
    LONGITUDE = ds.createVariable(varname='LONGITUDE', datatype='float', dimensions=("INSTRUMENT"))
    NOMINAL_DEPTH = ds.createVariable(varname='NOMINAL_DEPTH', datatype='float', dimensions=('INSTRUMENT'))

    ## main loop
    for index, file in enumerate(files_to_agg):
        print(index)
        with xr.open_dataset(file) as nc:

            ## in water data only
            nc = in_water(nc)

            start = sum(varlen_list[:index + 1])
            end = sum(varlen_list[:index + 2])
            n_cells = get_nvalues(nc)[1]
            UCUR[start:end] = get_varvalues(nc, 'UCUR')
            UCURqc[start:end] = get_varvalues(nc, 'UCUR_quality_control')
            VCUR[start:end] = get_varvalues(nc, 'VCUR')
            VCURqc[start:end] = get_varvalues(nc, 'VCUR_quality_control')
            if 'WCUR' in nc.data_vars:
                WCUR[start:end] = get_varvalues(nc, 'WCUR')
                WCURqc[start:end] = get_varvalues(nc, 'WCUR_quality_control')
            else:
                WCUR[start:end] = np.full(varlen_list[index], np.nan)
                WCURqc[start:end] = np.full(varlen_list[index], np.nan)
            ##calculate depth
            if 'HEIGHT_ABOVE_SENSOR' in nc.dims:
                DEPTH[start:end] = (nc.DEPTH + nc.HEIGHT_ABOVE_SENSOR.T).values.flatten()
                DEPTHqc[start:end] = np.array(n_cells * [nc.DEPTH_quality_control.values]).flatten()
            else:
                DEPTH[start:end] = nc.DEPTH.values
                DEPTHqc[start:end] = nc.DEPTH_quality_control.values
            ## set TIME and instrument index
            TIME[start:end] = (np.repeat(get_varvalues(nc, 'TIME'), n_cells) - epoch) / one_day
            instrument_index[start:end] = np.repeat(index, varlen_list[index + 1])
            ## get and store deployment metadata
            LATITUDE[index] = nc.LATITUDE.values
            LONGITUDE[index] = nc.LONGITUDE.values
            NOMINAL_DEPTH[index] = TStools.get_nominal_depth(nc)
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

    contributor_name, contributor_email, contributor_role = TStools.get_contributors(files_to_agg)
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
                    'contributor_role':        "; ".join(contributor_role)}

    global_attribute_dictionary.update(add_attribute)
    ds.setncatts(dict(sorted(global_attribute_dictionary.items())))

    ds.close()


    ## create the output file name and rename the tmp file
    facility_code = TStools.get_facility_code(files_to_agg[0])
    data_code = 'VZ'
    product_type = 'aggregated-timeseries'
    file_version = 1
    output_name = '_'.join(['IMOS', facility_code, data_code, time_start_filename, site_code, ('FV0'+str(file_version)),
                            ("velocity-"+product_type),
                            ('END-'+ time_end_filename), 'C-' + datetime.utcnow().strftime(file_timeformat)]) + '.nc'
    ncout_path = os.path.join(base_path, output_name)
    os.rename(outfile, ncout_path)

    return ncout_path, bad_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenate X,Y,Z velocity variables from ALL instruments from ALL deployments from ONE site")
    parser.add_argument('-site', dest='site_code', help='site code, like NRMMAI',  required=True)
    parser.add_argument('-files', dest='filenames', help='name of the file that contains the source URLs', required=True)
    parser.add_argument('-path', dest='output_path', help='path where the result file will be written. Default: ./', default='./', required=False)
    args = parser.parse_args()

    with open(args.filenames) as ff:
        files_to_agg = [line.rstrip() for line in ff]


    print(aggregate_velocity(files_to_agg=files_to_agg, site_code=args.site_code, base_path = args.output_path))
