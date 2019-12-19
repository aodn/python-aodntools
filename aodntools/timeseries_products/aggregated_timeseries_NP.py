from netCDF4 import Dataset
from dateutil.parser import parse
import numpy as np
import bisect
import argparse
import sys
import os.path
import json
from datetime import datetime
import glob
from pkg_resources import resource_filename

import xarray as xr
import pandas as pd

from aodntools import __version__


TEMPLATE_JSON = resource_filename(__name__, 'aggregated_timeseries_template.json')


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


def get_contributors(files_to_agg):
    """
    get the author and principal investigator details for each file

    :param files_to_aggregate: list of files
    :return: list: contributor_name, email and role
    """

    contributors = set()
    contributor_name, contributor_email, contributor_role = [], [], []

    for file in files_to_agg:
        with xr.open_dataset(file) as nc:
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

    agg_attr = {'title': ("Long Timeseries Aggregated product: " + varname + " at " + site + " between " + \
                          pd.to_datetime(agg_dataset.TIME.values.min()).strftime(timeformat) + " and " + \
                          pd.to_datetime(agg_dataset.TIME.values.max()).strftime(timeformat)),
                'site_code': site,
                'local_time_zone': '',
                'time_coverage_start': pd.to_datetime(agg_dataset.TIME.values.min()).strftime(timeformat),
                'time_coverage_end': pd.to_datetime(agg_dataset.TIME.values.max()).strftime(timeformat),
                'geospatial_vertical_min': float(agg_dataset.DEPTH.min()),
                'geospatial_vertical_max': float(agg_dataset.DEPTH.max()),
                'geospatial_lat_min': agg_dataset.LATITUDE.values.min(),
                'geospatial_lat_max': agg_dataset.LATITUDE.values.max(),
                'geospatial_lon_min': agg_dataset.LONGITUDE.values.min(),
                'geospatial_lon_max': agg_dataset.LONGITUDE.values.max(),
                'date_created': datetime.utcnow().strftime(timeformat),
                'history': datetime.utcnow().strftime(timeformat) + ': Aggregated file created.',
                'keywords': ', '.join(list(agg_dataset.variables) + ['AGGREGATED'])}
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
    if len(add_variable_attribute) > 0:
        for key in add_variable_attribute.keys():
            variable_attributes[key].update(add_variable_attribute[key])

    return variable_attributes


def get_data_code(VoI):
    """
    get data code sensu IMOS conventions from variable code

    :param VoI: variable code
    :return: variable data code
    """

    # dictionary of data code. could be read from external file
    dataCodes = {'DEPTH': 'Z',
                 'PRES': 'Z',
                 'PRES_REL': 'Z',
                 'TEMP': 'T',
                 'PSAL': 'S',
                 'PAR': 'F',
                 'TURB': 'U',
                 'TURBF': 'U',
                 'DOX1': 'O',
                 'DOX1_2': 'O',
                 'DOX1_3': 'O',
                 'DOX2': 'O',
                 'DOX2_1': 'O',
                 'DOXS': 'O',
                 'CPHL': 'B',
                 'CHLU': 'B',
                 'CHLF': 'B'}
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

    output_name = '_'.join(
        ['IMOS', facility_code, data_code, t_start, site_code, ('FV0' + str(file_version)), (VoI + "-" + product_type),
         ('END-' + t_end), 'C-' + datetime.utcnow().strftime(file_timeformat)]) + '.nc'

    return output_name


def create_empty_dataframe(columns):
    """
    create empty dataframe from a dict with data types

    :param: variable name and variable file. List of tuples
    :return: empty dataframe
    """

    return pd.DataFrame({k: pd.Series(dtype=t) for k, t in columns})


def write_netCDF_aggfile(agg_dataset, output_path, encoding):
    """
    write netcdf file

    :param agg_dataset: aggregated xarray dataset
    :param output_path: full path of the netCDF file to be written
    :return: name of the netCDf file written
    """

    agg_dataset.to_netcdf(output_path, encoding=encoding, format='NETCDF4_CLASSIC')

    return output_path


def in_water(nc):
    """
    cut data to in-water only timestamps, dropping resulting NaN.

    :param nc: xarray dataset
    :return: xarray dataset
    """
    time_deployment_start = np.datetime64(nc.attrs['time_deployment_start'][:-1])
    time_deployment_end = np.datetime64(nc.attrs['time_deployment_end'][:-1])
    TIME = nc['TIME'][:]
    return nc.where((TIME >= time_deployment_start) & (TIME <= time_deployment_end), drop=True)





#####################################################################################################
## MAIN FUNCTION
#####################################################################################################

def main_aggregator(files_to_aggregate, var_to_aggregate, site_code, base_path='./'):
    """
    Aggregates the variable of interest, its coordinates, quality control and metadata variables, from each file in
    the list into a netCDF file and returns its file name.

    :param files_to_agg: List of URLs for files to aggregate.
    :param var_to_agg: Name of variable to aggregate.
    :param site_code: code of the mooring site.
    :param base_path: path where the result file will be written
    :return: File path of the aggregated product
    :rtype: string
    """

    VoI = var_to_aggregate

    with open(TEMPLATE_JSON) as json_file:
        variable_attribute_dictionary = json.load(json_file)['_variables']


    ## sort the file URL in chronological order of deployment
    files_to_aggregate = sort_files_to_aggregate(files_to_aggregate)

    ## empty container for metadata
    metadata = pd.DataFrame(columns=['source_file', 'instrument_id', 'LATITUDE', 'LONGITUDE', 'NOMINAL_DEPTH'])

    ## empty np arrays
    time_all = np.array([], dtype='datetime64[ns]')
    VoI_all = np.array([], dtype='float64')
    depth_all = np.array([], dtype='float64')
    pressure_all = np.array([], dtype='float64')
    pressurerel_all = np.array([], dtype='float64')

    VoIQC_all = np.array([], dtype='int')
    depthQC_all = np.array([], dtype='int')
    pressureQC_all = np.array([], dtype='int')
    pressurerelQC_all = np.array([], dtype='int')

    instrument_index_all = np.array([], dtype='int')

    applied_offset = []  ## to store the PRES_REL attribute which could vary by deployment
    rejected_files = []
    bad_files = {}

    #print("Processing %i files..." % len(files_to_aggregate))

    for index, file in enumerate(files_to_aggregate):
        #print(index)
        try:
            with xr.open_dataset(file) as nc:
                file_problems = check_file(nc, VoI, site_code, variable_attribute_dictionary)
                if file_problems == []:
                    nc = in_water(nc)
                    ## get variable list
                    varList = list(nc.data_vars)

                    ## get metadata
                    metadata = metadata.append({'source_file': file,
                                                'instrument_id': nc.attrs['deployment_code'] + '; ' + nc.attrs[
                                                    'instrument'] + '; ' + nc.attrs['instrument_serial_number'],
                                                'LATITUDE': nc.LATITUDE.squeeze().values,
                                                'LONGITUDE': nc.LONGITUDE.squeeze().values,
                                                'NOMINAL_DEPTH': get_nominal_depth(nc)},
                                               ignore_index=True)

                    ## get the variable values
                    VoI_tmp = nc[VoI].values
                    VoIqc_tmp = nc[VoI + '_quality_control'].values

                    if 'DEPTH' in varList:
                        depth = nc.DEPTH.values
                        depthqc = nc.DEPTH_quality_control.values
                    else:
                        depth = np.full(len(VoI_tmp), np.nan)
                        depthqc = np.full(len(VoI_tmp), np.nan)

                    if 'PRES' in varList:
                        pressure = nc.PRES.values
                        pressureQC = nc.PRES_quality_control.values
                    else:
                        pressure = np.full(len(VoI_tmp), np.nan)
                        pressureQC = np.full(len(VoI_tmp), np.nan)

                    if 'PRES_REL' in varList:
                        pressureRel = nc.PRES_REL.values
                        pressureRelQC = nc.PRES_REL_quality_control.values
                        try:
                            applied_offset.append(nc.PRES_REL.applied_offset)
                        except:
                            applied_offset.append(np.nan)
                    else:
                        pressureRel = np.full(len(VoI_tmp), np.nan)
                        pressureRelQC = np.full(len(VoI_tmp), np.nan)
                        applied_offset.append(np.nan)

                    ## TIME and Instrument index
                    time = nc.TIME.values
                    instrument_index = np.array([index] * len(time))

                    ## Concatenate arrays
                    VoI_all = np.concatenate((VoI_all, VoI_tmp))
                    VoIQC_all = np.concatenate((VoIQC_all, VoIqc_tmp))

                    pressure_all = np.concatenate((pressure_all, pressure))
                    pressureQC_all = np.concatenate((pressureQC_all, pressureQC))
                    pressurerel_all = np.concatenate((pressurerel_all, pressureRel))
                    pressurerelQC_all = np.concatenate((pressurerelQC_all, pressureRelQC))

                    depth_all = np.concatenate((depth_all, depth))
                    depthQC_all = np.concatenate((depthQC_all, depthqc))

                    time_all = np.concatenate((time_all, time))
                    instrument_index_all = np.concatenate((instrument_index_all, instrument_index))
                else:
                    rejected_files.append(file)
                    bad_files.update({file: file_problems})
        except:
            print("UNKNOWN ERROR: " + file, file=sys.stderr)

    ds = xr.Dataset({VoI: (['OBSERVATION'], VoI_all),
                     VoI + '_quality_control': (['OBSERVATION'], VoIQC_all),
                     'DEPTH': (['OBSERVATION'], depth_all),
                     'DEPTH_quality_control': (['OBSERVATION'], depthQC_all),
                     'PRES': (['OBSERVATION'], pressure_all),
                     'PRES_quality_control': (['OBSERVATION'], pressureQC_all),
                     'PRES_REL': (['OBSERVATION'], pressurerel_all),
                     'PRES_REL_quality_control': (['OBSERVATION'], pressurerelQC_all),
                     'TIME': (['OBSERVATION'], time_all),
                     'instrument_index': (['OBSERVATION'], instrument_index_all)})

    metadata.index.rename('INSTRUMENT', inplace=True)

    ds = xr.merge([ds, xr.Dataset.from_dataframe(metadata)])
    ds = ds.drop('INSTRUMENT')      ## this is needed as the merge creates the dataframe index as a variable


    ## set variable attributes
    varlist = list(ds.data_vars)
    add_variable_attribute = {'PRES_REL': {'applied_offset_by_instrument': applied_offset}}
    variable_attributes = set_variableattr(varlist, variable_attribute_dictionary, add_variable_attribute)
    time_units = variable_attributes['TIME'].pop('units')
    time_calendar = variable_attributes['TIME'].pop('calendar')
    for variable in varlist:
        ds[variable].attrs = variable_attributes[variable]

    ## set global attrs
    contributor_name, contributor_email, contributor_role = get_contributors(files_to_aggregate)
    add_attribute = {'rejected_files': "\n".join(rejected_files),
                     'contributor_name': "; ".join(contributor_name),
                     'contributor_email': "; ".join(contributor_email),
                     'contributor_role': "; ".join(contributor_role)}
    ds.attrs = set_globalattr(ds, TEMPLATE_JSON, VoI, site_code, add_attribute)

    ## set compression & encoding
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    encoding.update({'TIME':                {'_FillValue': False,
                                             'units': time_units,
                                             'calendar': time_calendar},
                'LONGITUDE':                {'_FillValue': False},
                'LATITUDE':                 {'_FillValue': False},
                'instrument_id':            {'dtype': '|S256'},
                'source_file':              {'dtype': '|S256'}})


    ## create the output file name and write the aggregated product as netCDF
    facility_code = get_facility_code(files_to_aggregate[0])
    data_code = get_data_code(VoI) + 'Z'
    product_type = 'aggregated-timeseries'
    file_version = 1
    ncout_filename = generate_netcdf_output_filename(nc=ds, facility_code=facility_code, data_code=data_code,
                                                     VoI=VoI, site_code=site_code, product_type=product_type,
                                                     file_version=file_version)
    ncout_path = os.path.join(base_path, ncout_filename)
    write_netCDF_aggfile(ds, ncout_path, encoding)

    return ncout_path, bad_files


## END




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenate ONE variable from ALL instruments from ALL deployments from ONE site")
    parser.add_argument('-var', dest='varname', help='name of the variable to concatenate. Like TEMP, PSAL', required=True)
    parser.add_argument('-site', dest='site_code', help='site code, like NRMMAI',  required=True)
    parser.add_argument('-files', dest='filenames', help='name of the file that contains the source URLs', required=True)
    parser.add_argument('-path', dest='output_path', help='path where the result file will be written. Defaul ./', default='./', required=False)
    args = parser.parse_args()

    files_to_aggregate = pd.read_csv(args.filenames, header=None)[0].tolist()

    print(main_aggregator(files_to_aggregate=files_to_aggregate, var_to_aggregate=args.varname, site_code=args.site_code, base_path = args.output_path))

