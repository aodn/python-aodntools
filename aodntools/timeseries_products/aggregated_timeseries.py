from __future__ import print_function
import sys
import os.path
from dateutil.parser import parse
from datetime import datetime
import json
from netCDF4 import Dataset
import argparse

import numpy as np
import xarray as xr
import pandas as pd

from geoserverCatalog import get_moorings_urls



def sort_files_to_aggregate(files_to_agg):
    """
    sort the list of files to aggregate by time_deployment start attribute

    :param files_to_agg: list of file URLs
    :return: list of file URLs
    """
    file_list_dataframe = pd.DataFrame(columns=["url", "deployment_date"])
    for file in files_to_agg:
        with Dataset(file) as nc:
            try:
                file_list_dataframe = file_list_dataframe.append({'url': file,
                                                                'deployment_date': parse(nc.getncattr('time_deployment_start'))},
                                                                ignore_index=True)
            except:
                print(file)

    file_list_dataframe = file_list_dataframe.sort_values(by='deployment_date')

    return list(file_list_dataframe['url'])


def check_file(nc, VoI, site_code, variable_attribute_dictionary):
    """
    Return True the file pass all the following tests:
    VoI is present
    TIME is present
    LATITUDE is present
    LONGITUDE is present
    NOMINAL_DEPTH is present as variable or attribute
    file_version is FV01
    Return False if at least one of the tests fail
    if LATITUDE is a dimension has length 1
    if LONGITUDE is a dimension has length 1

    :param file: name of the netcdf file
    :param nc: xarray dataset
    :param VoI: string. Variable of Interest
    :param site_code: code of the mooring site
    :param variable_attribute_dictionary: dictionary of variable attributes
    :return: dictionary with the file name and list of failed tests
    """

    attributes = list(nc.attrs)
    variables = list(nc.variables)
    dimensions = list(nc.dims)
    allowed_dimensions = ['TIME', 'LATITUDE', 'LONGITUDE']
    error_list = []

    nc_site_code = getattr(nc, 'site_code', '')
    if nc_site_code != site_code:
        error_list.append('Wrong site_code: ' + nc_site_code)

    nc_file_version = getattr(nc, 'file_version', '')
    if 'Level 1' not in nc_file_version:
        error_list.append('Wrong file version: ' + nc_file_version)

    if 'TIME' not in variables:
        error_list.append('TIME variable missing')

    if 'LATITUDE' not in variables:
        error_list.append('LATITUDE variable missing')

    if 'LONGITUDE' not in variables:
        error_list.append('LONGITUDE variable missing')

    if 'NOMINAL_DEPTH' not in variables and 'instrument_nominal_depth' not in attributes:
        error_list.append('no NOMINAL_DEPTH')

    if VoI not in variables:
        error_list.append(VoI + ' not in file')
    else:
        VoIdimensions = list(nc[VoI].dims)
        if 'TIME' not in VoIdimensions:
            error_list.append('TIME is not a dimension')
        if 'LATITUDE' in VoIdimensions and len(nc.LATITUDE) > 1:
            error_list.append('more than one LATITUDE')
        if 'LONGITUDE' in VoIdimensions and len(nc.LONGITUDE) > 1:
            error_list.append('more than one LONGITUDE')
        for d in range(len(VoIdimensions)):
            if VoIdimensions[d] not in allowed_dimensions:
                error_list.append('not allowed dimensions: ' + VoIdimensions[d])
                break
        # test the units. not used. To be implemented in the future
        # if nc[VoI].attrs['units'] != variable_attribute_dictionary[VoI]['units']:
        #     error_list.append('Wrong units: ' + nc[VoI].attrs['units'])

    return error_list

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


def set_globalattr(agg_dataset, templatefile, varname, site, add_attribute):
    """
    global attributes from a reference nc file and nc file

    :param agg_dataset: aggregated xarray dataset
    :param templatefile: name of the attributes JSON file
    :param varname: name of the variable of interest to aggregate
    :param site: site code
    :param add_attribute: dictionary of additional attributes to add name:value
    :return: dictionary of global attributes
    """

    timeformat = '%Y-%m-%dT%H:%M:%SZ'
    with open(templatefile) as json_file:
        global_metadata = json.load(json_file)["_global"]

    agg_attr = {'title':                    ("Long Timeseries Aggregated product: " + varname + " at " + site + " between " + \
                                             pd.to_datetime(agg_dataset.TIME.values.min()).strftime(timeformat) + " and " + \
                                             pd.to_datetime(agg_dataset.TIME.values.max()).strftime(timeformat)),
                'site_code':                site,
                'local_time_zone':          '',
                'time_coverage_start':      pd.to_datetime(agg_dataset.TIME.values.min()).strftime(timeformat),
                'time_coverage_end':        pd.to_datetime(agg_dataset.TIME.values.max()).strftime(timeformat),
                'geospatial_vertical_min':  float(agg_dataset.DEPTH.min()),
                'geospatial_vertical_max':  float(agg_dataset.DEPTH.max()),
                'geospatial_lat_min':       agg_dataset.LATITUDE.values.min(),
                'geospatial_lat_max':       agg_dataset.LATITUDE.values.max(),
                'geospatial_lon_min':       agg_dataset.LONGITUDE.values.min(),
                'geospatial_lon_max':       agg_dataset.LONGITUDE.values.max(),
                'date_created':             datetime.utcnow().strftime(timeformat),
                'history':                  datetime.utcnow().strftime(timeformat) + ': Aggregated file created.',
                'keywords':                 ', '.join(list(agg_dataset.variables) + ['AGGREGATED'])}
    global_metadata.update(agg_attr)
    global_metadata.update(add_attribute)

    return dict(sorted(global_metadata.items()))


def set_variableattr(varlist, variable_attribute_dictionary, add_variable_attribute):
    """
    set variables variables atributes

    :param varlist: list of variable names
    :param variable_attribute_dictionary: dictionary of the variable attributes
    :param add_variable_attribute: additional attributes to add
    :return: dictionary of attributes
    """

    # with open(templatefile) as json_file:
    #     variable_metadata = json.load(json_file)['_variables']
    variable_attributes = {key: variable_attribute_dictionary[key] for key in varlist}
    if len(add_variable_attribute)>0:
        for key in add_variable_attribute.keys():
            variable_attributes[key].update(add_variable_attribute[key])

    return variable_attributes


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
                 'DOX1':        'O',
                 'DOX1_2':      'O',
                 'DOX1_3':      'O',
                 'DOX2':        'O',
                 'DOX2_1':      'O',
                 'DOXS':        'O',
                 'CPHL':        'B',
                 'CHLU':        'B',
                 'CHLF':        'B'}
    return dataCodes[VoI]

def get_facility_code(fileURL):
    """
    get the facility code from the file URL

    :param fileURL: URL of a file
    :return: facility code
    """

    return os.path.basename(fileURL).split("_")[1]


def generate_netcdf_output_filename(nc, facility_code, data_code, VoI, site_code, product_type, file_version):
    """
    generate the output filename for the VoI netCDF file

    :param nc: aggregated dataset
    :param facility_code: facility code from file name
    :param data_code: data code sensu IMOS convention
    :param VoI: name of the variable to aggregate
    :param product_type: name of the product
    :param file_version: version of the output file
    :return: name of the output file
    """

    file_timeformat = '%Y%m%d'

    if '_' in VoI:
        VoI = VoI.replace('_', '-')
    t_start = pd.to_datetime(nc.TIME.min().values).strftime(file_timeformat)
    t_end = pd.to_datetime(nc.TIME.max().values).strftime(file_timeformat)

    output_name = '_'.join(['IMOS', facility_code, data_code, t_start, site_code, ('FV0'+str(file_version)), (VoI+"-"+product_type), ('END-'+ t_end), 'C-' + datetime.utcnow().strftime(file_timeformat)]) + '.nc'

    return output_name

def create_empty_dataframe(columns):
    """
    create empty dataframe from a dict with data types

    :param: variable name and variable file. List of tuples
    :return: empty dataframe
    """

    return pd.DataFrame({k: pd.Series(dtype=t) for k, t in columns})


def write_netCDF_aggfile(agg_dataset, ncout_filename, encoding, base_path):
    """
    write netcdf file

    :param agg_dataset: aggregated xarray dataset
    :param ncout_filename: name of the netCDF file to be written
    :return: name of the netCDf file written
    """

    agg_dataset.to_netcdf(os.path.join(base_path, ncout_filename), encoding=encoding, format='NETCDF4_CLASSIC')

    return ncout_filename

## MAIN FUNCTION
def main_aggregator(files_to_agg, var_to_agg, site_code, base_path='./'):
    """
    Aggregates the variable of interest, its coordinates, quality control and metadata variables, from each file in
    the list into a netCDF file and returns its file name.

    :param files_to_agg: List of URLs for files to aggregate.
    :param var_to_agg: Name of variable to aggregate.
    :param site_code: code of the mooring site.
    :param base_path: path where the result file will be written
    :return: File name of the aggregated product
    :rtype: string
    """

    ## constants
    FILLVALUE = 999999.0
    variable_attributes_templatefile = 'TSagg_metadata.json'

    ## sort the file URL in chronological order of deployment
    files_to_agg = sort_files_to_aggregate(files_to_agg)

    var_to_agg_qc = var_to_agg + '_quality_control'
    ## create empty DF for main and auxiliary variables
    MainDF_types = [(var_to_agg, float),
                    (var_to_agg_qc, np.byte),
                    ('TIME', np.float64),
                    ('DEPTH', float),
                    ('DEPTH_quality_control', np.byte),
                    ('PRES', np.float64),
                    ('PRES_quality_control', np.byte),
                    ('PRES_REL', np.float64),
                    ('PRES_REL_quality_control', np.byte),
                    ('instrument_index', int)]

    AuxDF_types = [('source_file', str),
                   ('instrument_id', str),
                   ('LONGITUDE', float),
                   ('LATITUDE', float),
                   ('NOMINAL_DEPTH', float)]

    variableMainDF = create_empty_dataframe(MainDF_types)
    variableAuxDF = create_empty_dataframe(AuxDF_types)

    ## main loop

    ## get the variables attribute dictionary
    with open(variable_attributes_templatefile) as json_file:
        variable_attribute_dictionary = json.load(json_file)['_variables']

    fileIndex = 0
    rejected_files = []
    bad_files = {}
    applied_offset =[]      ## to store the PRES_REL attribute which could vary by deployment
    for file in files_to_agg:
        print(fileIndex, end=" ")
        sys.stdout.flush()

        try:
            ## it will open the netCDF files as a xarray Dataset
            with xr.open_dataset(file, decode_times=True) as nc:
                ## do only if the file pass all the sanity tests
                file_problems = check_file(nc, var_to_agg, site_code, variable_attribute_dictionary)
                if file_problems == []:
                    varnames = list(nc.variables.keys())
                    nobs = len(nc.TIME)

                    ## get the in-water times
                    ## important to remove the timezone aware of the converted datetime object from a string
                    time_deployment_start = pd.to_datetime(parse(nc.attrs['time_deployment_start'])).tz_localize(None)
                    time_deployment_end = pd.to_datetime(parse(nc.attrs['time_deployment_end'])).tz_localize(None)

                    DF = pd.DataFrame({ var_to_agg: nc[var_to_agg].squeeze(),
                                        var_to_agg_qc: nc[var_to_agg + '_quality_control'].squeeze(),
                                        'TIME': nc.TIME.squeeze(),
                                        'instrument_index': np.repeat(fileIndex, nobs)})

                    ## check for DEPTH/PRES variables in the nc and its qc flags
                    if 'DEPTH' in varnames:
                        DF['DEPTH'] = nc.DEPTH.squeeze()
                        if 'DEPTH_quality_control' in varnames:
                            DF['DEPTH_quality_control'] = nc.DEPTH_quality_control.squeeze()
                        else:
                            DF['DEPTH_quality_control'] = np.repeat(0, nobs)
                    else:
                        DF['DEPTH'] = np.repeat(FILLVALUE, nobs)
                        DF['DEPTH_quality_control'] = np.repeat(9, nobs)

                    if 'PRES' in varnames:
                        DF['PRES'] = nc.PRES.squeeze()
                        if 'PRES_quality_control' in varnames:
                            DF['PRES_quality_control'] = nc.PRES_quality_control.squeeze()
                        else:
                            DF['PRES_quality_control'] = np.repeat(0, nobs)
                    else:
                        DF['PRES'] = np.repeat(FILLVALUE, nobs)
                        DF['PRES_quality_control'] = np.repeat(9, nobs)

                    if 'PRES_REL' in varnames:
                        DF['PRES_REL'] = nc.PRES_REL.squeeze()
                        try:
                            applied_offset.append(nc.PRES_REL.applied_offset)
                        except:
                            applied_offset.append(np.nan)
                        if 'PRES_REL_quality_control' in varnames:
                            DF['PRES_REL_quality_control'] = nc.PRES_REL_quality_control.squeeze()
                        else:
                            DF['PRES_REL_quality_control'] = np.repeat(0, nobs)
                    else:
                        DF['PRES_REL'] = np.repeat(FILLVALUE, nobs)
                        DF['PRES_REL_quality_control'] = np.repeat(9, nobs)
                        applied_offset.append(np.nan)


                    ## select only in water data
                    DF = DF[(DF['TIME']>=time_deployment_start) & (DF['TIME']<=time_deployment_end)]

                    ## append data
                    variableMainDF = pd.concat([variableMainDF, DF], ignore_index=True, sort=False)


                    # append auxiliary data
                    variableAuxDF = variableAuxDF.append({'source_file': file,
                                                          'instrument_id': nc.attrs['deployment_code'] + '; ' + nc.attrs['instrument'] + '; ' + nc.attrs['instrument_serial_number'],
                                                          'LONGITUDE': nc.LONGITUDE.squeeze().values,
                                                          'LATITUDE': nc.LATITUDE.squeeze().values,
                                                          'NOMINAL_DEPTH': get_nominal_depth(nc)}, ignore_index = True)
                    fileIndex += 1
                else:
                    rejected_files.append(file)
                    bad_files.update({file: file_problems})
        except:
            print("FILE NOT FOUND: " + file, file=sys.stderr)

    print()


    ## rename indices
    variableAuxDF.index.rename('INSTRUMENT', inplace=True)
    variableMainDF.index.rename('OBSERVATION', inplace=True)

    ## get the list of variables
    varlist = list(variableMainDF.columns) + list(variableAuxDF.columns)


    ## set variable attributes
    add_variable_attribute = {'PRES_REL': {'applied_offset_by_instrument': applied_offset}}
    variable_attributes = set_variableattr(varlist, variable_attribute_dictionary, add_variable_attribute)
    time_units = variable_attributes['TIME'].pop('units')
    time_calendar = variable_attributes['TIME'].pop('calendar')

    ## build the output file
    agg_dataset = xr.Dataset({var_to_agg:                   (['OBSERVATION'],variableMainDF[var_to_agg].astype('float32'), variable_attributes[var_to_agg]),
                          var_to_agg + '_quality_control':  (['OBSERVATION'],variableMainDF[var_to_agg_qc].astype(np.byte), variable_attributes[var_to_agg+'_quality_control']),
                          'TIME':                           (['OBSERVATION'],variableMainDF['TIME'], variable_attributes['TIME']),
                          'DEPTH':                          (['OBSERVATION'],variableMainDF['DEPTH'].astype('float32'), variable_attributes['DEPTH']),
                          'DEPTH_quality_control':          (['OBSERVATION'],variableMainDF['DEPTH_quality_control'].astype(np.byte), variable_attributes['DEPTH_quality_control']),
                          'PRES':                           (['OBSERVATION'],variableMainDF['PRES'].astype('float32'), variable_attributes['PRES']),
                          'PRES_quality_control':           (['OBSERVATION'],variableMainDF['PRES_quality_control'].astype(np.byte), variable_attributes['PRES_quality_control']),
                          'PRES_REL':                       (['OBSERVATION'],variableMainDF['PRES_REL'].astype('float32'), variable_attributes['PRES_REL']),
                          'PRES_REL_quality_control':       (['OBSERVATION'],variableMainDF['PRES_REL_quality_control'].astype(np.byte), variable_attributes['PRES_REL_quality_control']),
                          'instrument_index':               (['OBSERVATION'],variableMainDF['instrument_index'].astype('int64'), variable_attributes['instrument_index']),
                          'LONGITUDE':                      (['INSTRUMENT'], variableAuxDF['LONGITUDE'].astype('float32'), variable_attributes['LONGITUDE']),
                          'LATITUDE':                       (['INSTRUMENT'], variableAuxDF['LATITUDE'].astype('float32'), variable_attributes['LATITUDE']),
                          'NOMINAL_DEPTH':                  (['INSTRUMENT'], variableAuxDF['NOMINAL_DEPTH']. astype('float32'), variable_attributes['NOMINAL_DEPTH']),
                          'instrument_id':                  (['INSTRUMENT'], variableAuxDF['instrument_id'].astype('|S256'), variable_attributes['instrument_id'] ),
                          'source_file':                    (['INSTRUMENT'], variableAuxDF['source_file'].astype('|S256'), variable_attributes['source_file'])})


    ## Set global attrs
    globalattr_file = 'TSagg_metadata.json'
    add_attribute = {'rejected_files': "\n".join(rejected_files)}
    agg_dataset.attrs = set_globalattr(agg_dataset, globalattr_file, var_to_agg, site_code, add_attribute)

    ## add version
    github_comment = ' Product created with https://github.com/aodn/data-services/blob/master/ANMN/LTSP/TSaggregator/aggregated_timeseries.py'

    agg_dataset.attrs['lineage'] += github_comment

    ## create the output file name and write the aggregated product as netCDF
    facility_code = get_facility_code(files_to_agg[0])
    data_code = get_data_code(var_to_agg) + 'Z'
    product_type='aggregated-time-series'
    file_version=1
    ncout_filename = generate_netcdf_output_filename(nc=agg_dataset, facility_code=facility_code, data_code=data_code, VoI=var_to_agg, site_code=site_code, product_type=product_type, file_version=file_version)

    encoding = {'TIME':                     {'_FillValue': False,
                                             'units': time_units,
                                             'calendar': time_calendar},
                'LONGITUDE':                {'_FillValue': False},
                'LATITUDE':                 {'_FillValue': False},
                'instrument_id':            {'dtype': '|S256'},
                'source_file':              {'dtype': '|S256'}}

    write_netCDF_aggfile(agg_dataset, ncout_filename, encoding, base_path)

    return ncout_filename, bad_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenate ONE variable from ALL instruments from ALL deployments from ONE site")
    parser.add_argument('-var', dest='varname', help='name of the variable to concatenate. Like TEMP, PSAL', required=True)
    parser.add_argument('-site', dest='site_code', help='site code, like NRMMAI',  required=True)
    parser.add_argument('-files', dest='filenames', help='name of the file that contains the source URLs', required=True)
    parser.add_argument('-path', dest='output_path', help='path where the result file will be written. Defaul ./', default='./', required=False)
    args = parser.parse_args()

    files_to_aggregate = pd.read_csv(args.filenames, header=None)[0].tolist()

    print(main_aggregator(files_to_agg=files_to_aggregate, var_to_agg=args.varname, site_code=args.site_code, base_path = args.output_path))
