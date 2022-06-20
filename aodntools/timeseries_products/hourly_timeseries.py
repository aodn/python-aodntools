#!/usr/bin/env python3

import argparse
import json
import os.path
from collections import OrderedDict

import numpy as np
import pandas as pd
import xarray as xr
from dateutil.parser import parse
from pkg_resources import resource_filename

from aodntools import __version__
from aodntools.timeseries_products import aggregated_timeseries as utils
from aodntools.timeseries_products.common import (NoInputFilesError, check_file, get_qc_variable_names, in_water,
                                                  current_utc_timestamp, TIMESTAMP_FORMAT, DATESTAMP_FORMAT)

TEMPLATE_JSON = resource_filename(__name__, 'hourly_timeseries_template.json')
BINNING_METHOD_JSON = resource_filename(__name__, 'binning_method.json')


def check_files(file_list, site_code, parameter_names_accepted, input_dir=''):
    """
    Return a chronologically sorted file_list and a dictionary if the file fails one or more of the tests

    :param file_list: list or file URLs
    :param site_code: code of the mooring site
    :param parameter_names_accepted: list of names of accepted parameters
    :param input_dir: base path where source files are stored
    :return: dictionary with the file name and list of failed tests, list good files chronologically ordered
    """

    file_list_dataframe = pd.DataFrame(columns=["url", "deployment_date"])
    error_dict = {}

    for file in file_list:
        with xr.open_dataset(os.path.join(input_dir, file)) as nc:
            error_list = check_file(nc, site_code, parameter_names_accepted)
            if error_list:
                error_dict.update({file: error_list})
            else:
                file_list_dataframe = file_list_dataframe.append({'url': file,
                                                                  'deployment_date': parse(nc.time_deployment_start)},
                                                                 ignore_index=True)

    file_list_dataframe = file_list_dataframe.sort_values(by='deployment_date')
    file_list = file_list_dataframe['url'].to_list()
    if file_list == []:
        raise NoInputFilesError("no valid input files to aggregate")

    return file_list, error_dict



def get_parameter_names(nc):
    """
    get the names of the parameters that HAVE _quality_control ancillary var
    remove from the list coordinates with QC flag variable

    :param nc: xarray dataset
    :return: list of names
    """
    params = list(set([s.strip('_quality_control') for s in get_qc_variable_names(nc)]) - set(list(nc.coords)))
    return params


def good_data_only(nc, qcflags):
    """
    mask all the variables with QC for QC value less or equal than the specified

    :param nc: xarray dataset
    :param qcflags: list of QCflags indicating what variables to keep
    :return: xarray masked Dataset, dictionary of % of qced values per variable
    """
    varnames = get_parameter_names(nc)
    nc_masked = nc[varnames[0]].where(nc[varnames[0] + '_quality_control'].isin(qcflags)).to_dataset(name=varnames[0])
    for variable in varnames[1:]:
        nc_masked[variable] = nc[variable].where(nc[variable + '_quality_control'].isin(qcflags))

    return nc_masked


def get_QCcount (nc, qcflags):
    """
    count the number of qced values in the file
    :param nc: xarray dataset
    :param qcflags: QCflags to count
    :return: dictionary with % of registers QCed
    """
    qc_total_count = {}
    if 0 in qcflags and len(qcflags)>1:
        varnames = get_parameter_names(nc)

        for variable in varnames:
            flag_count = []
            for flag in qcflags:
                flag_count.append(int(np.sum(nc[variable+'_quality_control']==flag)))
            qc_total_count[variable] = {'qc0_count': flag_count[0]}
            qc_total_count[variable].update({'qcnon0_count': sum(flag_count[1:])})

    return qc_total_count

def update_QCcount(qc_count_all, qc_count):
    """
    Update qc count dictionary
    :param qc_count_all: dictionary with all variables qc count to be updated
    :param qc_count: dictionary of all variable qc count from one file
    :return: dictionary of qc proportions per variable
    """
    for variable in qc_count.keys():
        if variable in qc_count_all.keys():
            qc_count_all[variable]['qc0_count'] += qc_count[variable]['qc0_count']
            qc_count_all[variable]['qcnon0_count'] += qc_count[variable]['qcnon0_count']
        else:
            qc_count_all[variable] = qc_count[variable]

    return qc_count_all

def get_QC_percent(qc_count):
    """
    Calculate the % of qc values in the variables of a file
    :param qc_count: dictionary of qc counts
    :return: dictionary of % of qc values per variable
    """
    qc_percent = {}
    if len(qc_count) > 0:
        for variable in qc_count.keys():
            if qc_count[variable]['qcnon0_count'] > 0:
                qc_percent[variable] = {'percent_quality_controlled': round(100*(1-qc_count[variable]['qc0_count']/(qc_count[variable]['qcnon0_count'] + qc_count[variable]['qc0_count'])),2)}
            else:
                qc_percent[variable] = {'percent_quality_controlled': 0.00}

    return qc_percent



def get_nominal_depth(nc):
    """
    return nominal depth from NOMINAL_DEPTH variable or
    if it is not present from instrument_nominal_depth global attribute

    :param nc: xarray dataset
    :return: nominal depth of the instrument
    """

    if 'NOMINAL_DEPTH' in list(nc.variables):
        nominal_depth = nc.NOMINAL_DEPTH.squeeze().values
    else:
        nominal_depth = nc.instrument_nominal_depth

    return nominal_depth


def set_globalattr(nc_aggregated, templatefile, site_code, add_attribute, parameter_names):
    """
    global attributes from a reference nc file and nc file

    :param nc_aggregated: aggregated xarray dataset
    :param templatefile: name of the attributes JSON file
    :param site_code: code of the mooring site
    :param add_attribute: dictionary of additional attributes to add name:value
    :param parameter_names: list of aggregated parameters
    :return: dictionary of global attributes
    """

    timeformat = '%Y-%m-%dT%H:%M:%SZ'
    with open(templatefile) as json_file:
        global_metadata = json.load(json_file)["_global"]

    agg_attr = {'title': ("Long time series Hourly Aggregated product: all available non-velocity variables at " +
                          site_code + " between " + pd.to_datetime(nc_aggregated.TIME.values.min()).strftime(timeformat) +
                          " and " + pd.to_datetime(nc_aggregated.TIME.values.max()).strftime(timeformat)),
                'site_code': site_code,
                'time_coverage_start': pd.to_datetime(nc_aggregated.TIME.values.min()).strftime(timeformat),
                'time_coverage_end': pd.to_datetime(nc_aggregated.TIME.values.max()).strftime(timeformat),
                'geospatial_vertical_min': float(nc_aggregated.DEPTH.min()),
                'geospatial_vertical_max': float(nc_aggregated.DEPTH.max()),
                'geospatial_lat_min': nc_aggregated.LATITUDE.values.min(),
                'geospatial_lat_max': nc_aggregated.LATITUDE.values.max(),
                'geospatial_lon_min': nc_aggregated.LONGITUDE.values.min(),
                'geospatial_lon_max': nc_aggregated.LONGITUDE.values.max(),
                'date_created': current_utc_timestamp(),
                'history': current_utc_timestamp() + ': Hourly aggregated file created.',
                'keywords': ', '.join(parameter_names + ['HOURLY', 'AGGREGATED'])}
    global_metadata.update(agg_attr)
    global_metadata.update(add_attribute)

    return OrderedDict(sorted(global_metadata.items()))


def get_data_code(VoI):
    """
    get data code sensu IMOS conventions from variable code

    :param VoI: variable code
    :return: variable data code
    """

    # dictionary of data code. better if it is read from external file
    dataCodes = {'DEPTH': 'Z',
                 'PRES': 'Z',
                 'PRES_REL': 'Z',
                 'TEMP': 'T',
                 'PSAL': 'S',
                 'CNDC': 'C',
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

    if VoI in dataCodes:
        dCode = dataCodes[VoI]
    else:
        dCode = ""
    return dCode


def get_facility_code(fileURL):
    """
    get the facility code from the file URL

    :param fileURL: URL of a file
    :return: facility code
    """

    return os.path.basename(fileURL).split("_")[1]


def create_empty_dataframe(columns):
    """
    create empty dataframe from a dict with data types

    :param: variable name and variable file. List of tuples
    :return: empty dataframe
    """

    return pd.DataFrame({k: pd.Series(dtype=t) for k, t in columns})


def generate_netcdf_output_filename(nc, facility_code, data_code, site_code, product_type, file_version):
    """
    generate the output filename for the VoI netCDF file

    :param nc: aggregated dataset
    :param facility_code: facility code from file name
    :param data_code: data code sensu IMOS convention
    :param site_code: site code
    :param product_type: name of the product
    :param file_version: version of the output file
    :return: name of the output file
    """

    t_start = pd.to_datetime(nc.TIME.min().values).strftime(DATESTAMP_FORMAT)
    t_end = pd.to_datetime(nc.TIME.max().values).strftime(DATESTAMP_FORMAT)

    output_name = '_'.join(
        ['IMOS', facility_code, data_code, t_start, site_code, ('FV0' + str(file_version)), product_type,
         ('END-' + t_end), 'C-' + current_utc_timestamp(DATESTAMP_FORMAT)]) + '.nc'

    return output_name



def write_netCDF_aggfile(nc_aggregated, ncout_filename, encoding):
    """
    write netcdf file

    :param nc_aggregated: aggregated xarray dataset
    :param ncout_filename: path/name of the netCDF file to be written
    :param encoding: encoding dictionary
    :return: name of the netCDf file written
    """
    ## sort the variables in the data set
    variables_all = list(nc_aggregated.variables)
    variables_head = ['instrument_index', 'instrument_id', 'source_file', 'TIME', 'LONGITUDE', 'LATITUDE',
                      'NOMINAL_DEPTH', 'DEPTH', 'DEPTH_count', 'DEPTH_min', 'DEPTH_max', 'DEPTH_std', ]
    variables_rest = sorted(list(set(variables_all) - set(variables_head)))
    variables_all = variables_head + variables_rest

    nc_aggregated[variables_all].to_netcdf(ncout_filename, encoding=encoding,
                                           format='NETCDF4_CLASSIC')
    return ncout_filename


def PDresample_by_hour(df, function_dict, function_stats):
    """
    resample a dataframe by hour and calculate aggregation statistics
    (mean, std, min, max, count)
    the variables are renamed adding the corresponding suffix of the stats calculated
    from the function_names list

    :param function_stats: dictionary of binning ancillary stat functions
    :param function_dict: dictionary of binning methods to be applied to each variable
    :param df: pandas dataframe with ancillary variables and coords removed but with TIME as index
    :return: pandas dataframe
    """

    varnames = df.columns
    df_data = pd.DataFrame(index=pd.DatetimeIndex([]))
    for variable in varnames:
        ds_var = df[variable]
        ds_var_resample = ds_var.resample('1H', base=0.5)  # shift by half hour to centre bin on the hour
        ds_var_mean = ds_var_resample.apply(function_dict[variable]).astype(np.float32)
        df_data = pd.concat([df_data, ds_var_mean], axis=1, sort=False)
        for stat_method in function_stats:
            ds_var_stat = ds_var_resample.apply(stat_method).astype(np.float32)
            ds_var_stat = ds_var_stat.rename("_".join([variable, stat_method]))
            df_data = pd.concat([df_data, ds_var_stat], axis=1, sort=False)

    ##forward the index 30min so the timestamps are on the hour
    df_data.index += pd.to_timedelta('30min')

    return df_data




### MAIN FUNCTION
def hourly_aggregator(files_to_aggregate, site_code, qcflags, input_dir='', output_dir='./',
                      download_url_prefix=None, opendap_url_prefix=None):
    """
    Aggregate a dataset into 1 hour intervals and calculate related statistics

    :param files_to_aggregate: List of files to aggregate. Each path is interpreted relative
    to input_dir (if specified). These paths are listed in the `source_files` variable in the output file.
    :param site_code: code of the mooring site
    :param qcflags: list of QCflags indicating what values of the variables to keep
    :param input_dir: base path where source files are stored
    :param output_dir: path where the result file will be written
    :param download_url_prefix: URL prefix for file download (to be prepended to paths in files_to_aggregate)
    :param opendap_url_prefix: URL prefix for OPENAP access (to be prepended to paths in files_to_aggregate)
    :return: tuple (path of the output file, dict of rejected files and error messages)
    """

    parameter_names_accepted = ['DEPTH', 'CPHL', 'CHLF', 'CHLU', 'DOX', 'DOX1', 'DOX1_2', 'DOX1_3', 'DOX2',
                                'DOX2_1', 'DOXS', 'DOXY', 'PRES', 'PRES_REL', 'PSAL', 'TEMP', 'TURB', 'TURBF', 'PAR']
    function_stats = ['min', 'max', 'std', 'count']
    qcflags_names = {0: 'No_QC_performed', 1: 'Good_data', 2: 'Probably_good_data',
                     3: 'Bad_data_that_are_potentially_correctable', 4: 'Bad_data'}



    ## make sure that the list of qflags is sorted
    qcflags = sorted(qcflags)

    # Check files and sort chronologically
    files_to_aggregate, bad_files = check_files(files_to_aggregate, site_code, parameter_names_accepted,
                                                input_dir=input_dir)

    ## get binning function dictionary
    with open(BINNING_METHOD_JSON) as json_file:
        function_dict = json.load(json_file)


    ## get the variables attribute dictionary
    with open(TEMPLATE_JSON) as json_file:
        variable_attribute_dictionary = json.load(json_file)['_variables']

    df_data = pd.DataFrame()


    ## create empty DF with dtypes
    metadata_df_types = [('source_file', str),
                         ('instrument_id', str),
                         ('LONGITUDE', float),
                         ('LATITUDE', float),
                         ('NOMINAL_DEPTH', float)]
    df_metadata = create_empty_dataframe(metadata_df_types)

    ## containers
    parameter_names_all = []
    applied_offset = []
    qc_count_all = {}

    for file_index, file in enumerate(files_to_aggregate):
        print(file_index)
        with xr.open_dataset(os.path.join(input_dir, file), mask_and_scale=True, decode_times=True) as nc:
            parameter_names = list(set(list(nc.variables)) & set(parameter_names_accepted))
            parameter_names_all += parameter_names

            ## get PRES_REl offset, if exits
            if 'PRES_REL' in parameter_names:
                if 'applied_offset' in nc.PRES_REL.attrs:
                    applied_offset.append(nc.PRES_REL.applied_offset)
                else:
                    applied_offset.append(np.nan)

            nc_clean = in_water(nc)  # in water only
            qc_count = get_QCcount(nc_clean, qcflags)
            qc_count_all = update_QCcount(qc_count_all, qc_count)
            nc_clean = good_data_only(nc_clean, qcflags)  # good quality data only
            df_metadata = df_metadata.append({'source_file': file,
                                              'instrument_id': utils.get_instrument_id(nc),
                                              'LONGITUDE': nc.LONGITUDE.squeeze().values,
                                              'LATITUDE': nc.LATITUDE.squeeze().values,
                                              'NOMINAL_DEPTH': get_nominal_depth(nc)},
                                             ignore_index=True)

            # If TIME had out-of-range values before cleaning, nc_clean would now have a CFTimeIndex, which
            # breaks the resampling further down. Here we reset it to a DatetimeIndex as suggested here:
            # https://stackoverflow.com/questions/55786995/converting-cftime-datetimejulian-to-datetime/55787899#55787899
            if isinstance(nc_clean.indexes['TIME'], xr.coding.cftimeindex.CFTimeIndex):
                nc_clean['TIME'] = nc_clean.indexes['TIME'].to_datetimeindex()
            df_temp = nc_clean[parameter_names].to_dataframe()

            ## keep TIME as the only index (for ADCP files it would be a MultiIndex at this point)
            df_temp.reset_index(inplace=True)
            df_temp.set_index('TIME', inplace=True)
            df_temp = df_temp[parameter_names]

            df_temp = PDresample_by_hour(df_temp, function_dict, function_stats)  # do the magic
            df_temp['instrument_index'] = np.repeat(file_index, len(df_temp)).astype(np.int32)
            df_data = pd.concat([df_data, df_temp.reset_index()], ignore_index=True, sort=False)

    df_metadata.index.rename('INSTRUMENT', inplace=True)
    df_data.index.rename('OBSERVATION', inplace=True)
    ## rename index to TIME
    df_data.rename(columns={'index': 'TIME'}, inplace=True)

    qc_proportion_all = get_QC_percent(qc_count_all)


    nc_metadata = xr.Dataset({'LONGITUDE': (['INSTRUMENT'], df_metadata['LONGITUDE'].astype('float64')),
                              'LATITUDE': (['INSTRUMENT'], df_metadata['LATITUDE'].astype('float64')),
                              'NOMINAL_DEPTH': (['INSTRUMENT'], df_metadata['NOMINAL_DEPTH'].astype('float32')),
                              'instrument_id': (['INSTRUMENT'], df_metadata['instrument_id'].astype('|S256')),
                              'source_file': (['INSTRUMENT'], df_metadata['source_file'].astype('|S256'))})

    ## Check and drop all nan columns
    column_remove_list = []
    parameter_remove_list = []
    for parameter in parameter_names_all:
        if df_data[parameter].isna().all():
            column_remove_list.append(parameter)
            parameter_remove_list.append(parameter)
            for method in function_stats:
                column_remove_list.append(parameter+'_'+method)
    df_data.drop(columns=column_remove_list, inplace=True)
    ## remove the drop names from the parameter_names_all list
    parameter_names_all = list(set(parameter_names_all) - set(parameter_remove_list))

    nc_data = xr.Dataset.from_dataframe(df_data)
    nc_aggregated = xr.merge([nc_metadata, nc_data])
    nc_aggregated = nc_aggregated.drop('OBSERVATION')

    ## add global attributes
    add_attribute = {'rejected_files': "\n".join(list(bad_files)),
                     'included_values_flagged_as':  ", ".join([qcflags_names[flag] for flag in qcflags]),
                     'generating_code_version': __version__
                     }
    add_attribute.update(utils.get_contributors(files_to_aggregate, input_dir=input_dir))

    nc_aggregated.attrs = set_globalattr(nc_aggregated, TEMPLATE_JSON, site_code, add_attribute, parameter_names)
    nc_aggregated.attrs['abstract'] = nc_aggregated.attrs['abstract'].format(
        site_code=site_code,
        flags=", ".join(qcflags_names[flag] for flag in qcflags)
    )
    if 0 in qcflags:
        nc_aggregated.attrs['lineage'] += ('The percentage of quality controlled values used in the aggregation is '
                                           'indicated in the percent_quality_controlled variable attribute.')
    github_comment = ('\nThis file was created using https://github.com/aodn/python-aodntools/blob/'
                      '{v}/aodntools/timeseries_products/hourly_timeseries.py'.format(v=__version__)
                      )
    nc_aggregated.attrs['lineage'] += github_comment

    ## add variable attributes
    variablenames_others = ['TIME', 'LONGITUDE', 'LATITUDE', 'NOMINAL_DEPTH',
                            'instrument_index', 'instrument_id', 'source_file']
    parameter_names_all = list(set(parameter_names_all))
    variable_attributes = variable_attribute_dictionary
    variable_attributes['PRES_REL'].update({'applied_offset_by_instrument': applied_offset})
    if download_url_prefix or opendap_url_prefix:
        variable_attributes['source_file'].update(utils.source_file_attributes(download_url_prefix, opendap_url_prefix))

    time_units = variable_attributes['TIME'].pop('units')
    time_calendar = variable_attributes['TIME'].pop('calendar')
    encoding = {'TIME': {'_FillValue': None,
                         'units': time_units,
                         'calendar': time_calendar},
                'LONGITUDE': {'_FillValue': None},
                'LATITUDE': {'_FillValue': None},
                'NOMINAL_DEPTH': {'_FillValue': None},
                'instrument_id': {'dtype': '|S256'},
                'source_file': {'dtype': '|S256'}}


## add attributes to TIME, LAT/LON, and index variables
    for variable in variablenames_others:
        nc_aggregated[variable].attrs = variable_attributes[variable]

    for variable in parameter_names_all:
        ancillary_variables_attr = []
        ## remove the _FillValue attribute
        fill_value = variable_attributes[variable].pop('_FillValue')
        encoding.update({variable: {'_FillValue': fill_value}})
        ## replace nan by FillValue
        nc_aggregated[variable] = nc_aggregated[variable].fillna(fill_value)

        nc_aggregated[variable].attrs = variable_attributes[variable]
        nc_aggregated[variable].attrs['long_name'] = function_dict[variable] + " " + nc_aggregated[variable].attrs['long_name']
        nc_aggregated[variable].attrs.update({'cell_methods': 'TIME:' + function_dict[variable] + ' (interval: 1 hr comment: time mid point)'})

        ## add percent of QCed values
        if qc_proportion_all:
            nc_aggregated[variable].attrs.update(qc_proportion_all[variable])

        for stat_method in function_stats:
            variable_stat_name = variable + "_" + stat_method
            ancillary_variables_attr += [variable_stat_name]
            if stat_method == 'count':
                nc_aggregated[variable_stat_name].attrs['units'] = "1"
            else:
                nc_aggregated[variable_stat_name].attrs['units'] = variable_attributes[variable]['units']

            if 'standard_name' in nc_aggregated[variable].attrs:
                nc_aggregated[variable_stat_name].attrs['standard_name'] = nc_aggregated[variable].attrs['standard_name']
                nc_aggregated[variable+'_count'].attrs['standard_name'] = nc_aggregated[variable].attrs['standard_name'] + ' number_of_observations'

            nc_aggregated[variable_stat_name].attrs['long_name'] = stat_method + ' data value in the bin, after rejection of flagged data'
            nc_aggregated[variable_stat_name].attrs['cell_methods'] = 'TIME:' +  stat_method
            nc_aggregated[variable_stat_name].attrs['_FillValue'] = fill_value
            nc_aggregated[variable_stat_name] = nc_aggregated[variable_stat_name].fillna(fill_value)
        nc_aggregated[variable].attrs.update({'ancillary_variables': " ".join(ancillary_variables_attr)})


    ## create the output file name and write the aggregated product as netCDF
    facility_code = get_facility_code(os.path.join(input_dir, files_to_aggregate[0]))
    data_code = "".join(sorted(set(get_data_code(p) for p in parameter_names_all)))
    if 0 in qcflags:
        product_type = 'hourly-timeseries-including-non-QC'
    else:
        product_type = 'hourly-timeseries'
    file_version = 2
    ncout_filename = generate_netcdf_output_filename(nc=nc_aggregated, facility_code=facility_code, data_code=data_code,
                                                     site_code=site_code,
                                                     product_type=product_type, file_version=file_version)
    ncout_path = os.path.join(output_dir, ncout_filename)
    write_netCDF_aggfile(nc_aggregated, ncout_path, encoding)


    return ncout_path, bad_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenate ALL variables from ALL instruments from ALL deployments from ONE site at 1hr time bin")
    parser.add_argument('-site', dest='site_code', help='site code, like NRMMAI',  required=True)
    parser.add_argument('-files', dest='filenames', help='name of the file that contains the source URLs', required=True)
    parser.add_argument('-qc', dest='qcflags', help='list of QC flags to select variable values to keep', nargs='+', required=True)
    parser.add_argument('-indir', dest='input_dir', help='base path of input files', default='', required=False)
    parser.add_argument('-outdir', dest='output_dir', help='path where the result file will be written. Default ./',
                        default='./', required=False)
    args = parser.parse_args()

    with open(args.filenames, 'r') as file:
        files_to_aggregate = [i.strip() for i in file.readlines()]
    qcflags = [int(i) for i in args.qcflags]

    hourly_aggregator(files_to_aggregate=files_to_aggregate, site_code=args.site_code, qcflags=qcflags,
                      input_dir=args.input_dir, output_dir=args.output_dir)
