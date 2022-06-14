import numpy as np
import bisect
import argparse
import os.path
import json
from datetime import datetime, timezone

import xarray as xr
import pandas as pd

from pkg_resources import resource_filename

from aodntools import __version__
from aodntools.timeseries_products.common import current_utc_timestamp, TIMESTAMP_FORMAT, DATESTAMP_FORMAT
import aodntools.timeseries_products.aggregated_timeseries as TStools


TEMPLATE_JSON = resource_filename(__name__, 'gridded_timeseries_template.json')


def make_depth_bins(nc, increment=10):
    """
    generate a list of depth bins from min/max DEPTH.
    Re-scale the min/max to nearest increment round number and start the bin from there
    :param nc: hourly dataset
    :return: list of depth bins
    """

    dmin = int(nc.DEPTH.min() / increment) * increment
    dmax = int(nc.DEPTH.max() / increment) * increment + increment
    dbin = list(range(dmin, dmax, increment))

    return dbin

def adjust_depth_bins(depth_bins, depth_min, depth_max):
    """
    Adjust the provided depth bins to match the min and max registered depth
    :param depth_bins: list of target depths
    :param depth_min: min valid depth
    :param depth_max: max valid depth
    :return: list of adjusted target depths
    """

    target_depth_min = bisect.bisect_left(depth_bins, depth_min)
    target_depth_max = bisect.bisect_right(depth_bins, depth_max)
    return depth_bins[target_depth_min:target_depth_max]


def get_depth_mask(depth_bins, depths, max_separation):
    """
    return a boolean mask of depth bins where the interpolation is possible
    based on maximum separation between the instruments
    :param depth_bins: list of depth bins where the interpolation will happen
    :param depths: list of depth of the instruments at the current timestamp
    :param max_separation: maximum allowable separation between readings
    :return: list of depth masked
    """
    depth_mask = [False] * len(depth_bins)
    for depth_index in range(len(depths) - 1):
        if depths[depth_index + 1] - depths[depth_index] < max_separation:
            depth_min = bisect.bisect_left(depth_bins, depths[depth_index])
            depth_max = bisect.bisect_right(depth_bins, depths[depth_index + 1])
            depth_interval_len = depth_max - depth_min
            depth_mask[depth_min:depth_max] = [True] * depth_interval_len
    return depth_mask



def sort_depths(depths, values):
    """
    Sort the list of depths and values
    """
    index =list(range(len(depths)))
    index.sort(key=depths.__getitem__)
    sorted_depths = [depths[i] for i in index]
    sorted_values = [values[i] for i in index]

    return sorted_depths, sorted_values

def write_netCDF_aggfile(agg_dataset, output_path, encoding):
    """
    write netcdf file

    :param agg_dataset: aggregated xarray dataset
    :param output_path: full path of the netCDF file to be written
    :return: name of the netCDf file written
    """

    agg_dataset.to_netcdf(output_path, encoding=encoding, format='NETCDF4_CLASSIC')

    return output_path


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

    if '_' in VoI:
        VoI = VoI.replace('_', '-')
    t_start = pd.to_datetime(nc.TIME.min().values).strftime(DATESTAMP_FORMAT)
    t_end = pd.to_datetime(nc.TIME.max().values).strftime(DATESTAMP_FORMAT)

    output_name = '_'.join(['IMOS', facility_code, data_code, t_start, site_code, ('FV0'+str(file_version)), (VoI+"-"+product_type), ('END-'+ t_end), 'C-' + current_utc_timestamp(DATESTAMP_FORMAT)]) + '.nc'

    return output_name



## MAIN FUNCTION
def grid_variable(input_file, VoI, depth_bins=None, max_separation=50, depth_bins_increment=10,
                  input_dir='', output_dir='.', download_url_prefix=None, opendap_url_prefix=None):
    """
    Grid VoI into depth_bins.
    :param input_file: Input hourly aggregated file with VoI, DEPTH and TIME only (path interpreted relative
    to input_dir, if specified)
    :param VoI: variable of interest (TEMP or PSAL)
    :param depth_bins: list of depth where to interpolate. if null list is provided it will be calculated from the data
    :param max_separation: max separation allowed for instruments
    :param depth_bins_increment: in case no depth bins provided this is the increment for the calculated bins
    :param input_dir: base path where source files are stored
    :param output_dir: path where the result file will be written
    :param download_url_prefix: URL prefix for file download (to be prepended to input_file path)
    :param opendap_url_prefix: URL prefix for OPENAP access (to be prepended to input_file path)
    :return: path of interpolated output file
    """

    with xr.open_dataset(os.path.join(input_dir, input_file)) as nc_full:
        nc = nc_full[[VoI, 'TIME', 'DEPTH']]
        ## get lat/lon
        longitude_mean = nc_full.LONGITUDE.mean()
        latitude_mean = nc_full.LATITUDE.mean()

    site_code = nc.site_code
    ## get global attributes
    input_global_attributes = nc.attrs



    ## in case no depth bins provided, create depth bins to the nearest rounded depth
    ## if provided, adjust to the min-max registered depth
    if not depth_bins:
        depth_bins = make_depth_bins(nc=nc, increment=depth_bins_increment)
    else:
        depth_bins = [float(depth) for depth in depth_bins]     # in case depth_bins provided through inline arguments
        depth_bins = adjust_depth_bins(depth_bins, nc.DEPTH.min(), nc.DEPTH.max())

    ## create empty containers
    time_0 = pd.to_datetime('1950-01-01T00:00:00')
    time_min = nc.TIME.values.min()
    depth_bin_len = len(depth_bins)

    ## create empty containers for the interpolated values
    VoI_temp = xr.DataArray(np.full((depth_bin_len, 1), np.nan, dtype=np.float32), coords=[depth_bins, [time_0]],
                              dims=['DEPTH', 'TIME'])
    VoI_ndepths = xr.DataArray(np.full(1, 0, dtype='int'), coords=[[time_0]], dims=['TIME'])

    ## group nc by individual timestamps
    VoI_grouped = nc.groupby('TIME')


    for timestamp, group in VoI_grouped:
        time = [timestamp]
        n_depths = np.array(len(group[VoI]), dtype='int')

        if n_depths >= 2:
            VoI_values = list(group[VoI].values)
            depth = list(group.DEPTH.values)
            ## sort depths
            depth, VoI_values = sort_depths(depth, VoI_values)

            ## check for max separation
            depth_mask = get_depth_mask(depth_bins=depth_bins, depths=depth, max_separation=max_separation)
            ## do the interpolation
            interpolated_var = np.interp(depth_bins, depth, VoI_values, left=np.nan, right=np.nan)
            ## set masked depth bins to zero
            interpolated_var = interpolated_var * depth_mask
            interpolated_var[interpolated_var == 0] = np.nan

        else:
            interpolated_var = np.full((depth_bin_len, 1), np.nan)

        VoI_temp_tmp = xr.DataArray(interpolated_var.reshape(depth_bin_len, 1), coords=[depth_bins, time],
                                dims=['DEPTH', 'TIME'])
        VoI_ndepths_tmp = xr.DataArray([n_depths], coords=[time], dims=['TIME'])

        ## concatenate the interpolated values
        VoI_temp = xr.concat([VoI_temp, VoI_temp_tmp], dim='TIME')
        VoI_ndepths = xr.concat([VoI_ndepths, VoI_ndepths_tmp], dim='TIME')

    VoI_interpolated = xr.Dataset({VoI: VoI_temp.astype(np.float32),
                                   VoI + '_count': VoI_ndepths.astype('int16')})

    ## drop the very first record as it is dummy
    VoI_interpolated = VoI_interpolated.where(VoI_interpolated.TIME >= time_min, drop=True)

    ## Add lat/lon as scalar variables
    VoI_interpolated = VoI_interpolated.assign(LONGITUDE = longitude_mean,
                                               LATITUDE = latitude_mean)

    ## transpose dimensions to make CF compliant
    VoI_interpolated = VoI_interpolated.transpose('TIME', 'DEPTH')

    ## get the variables attribute dictionary
    with open(TEMPLATE_JSON) as json_file:
        attr_dictionary = json.load(json_file)

    variable_attribute_dictionary = attr_dictionary['_variables']
    global_attribute_dictionary = attr_dictionary['_global']

    ## set variable attributes
    varlist = list(VoI_interpolated.variables)
    add_variable_attribute = {}
    variable_attributes = set_variableattr(varlist, variable_attribute_dictionary, add_variable_attribute)
    time_units = variable_attributes['TIME'].pop('units')
    time_calendar = variable_attributes['TIME'].pop('calendar')
    for variable in varlist:
        VoI_interpolated[variable].attrs = variable_attributes[variable]

    ## set global attributes
    # copy selected attributes from input file
    for attr in ('geospatial_lat_min', 'geospatial_lat_max', 'geospatial_lon_min', 'geospatial_lon_max', 'site_code',
                 'included_values_flagged_as', 'contributor_name', 'contributor_role', 'contributor_email'):
        VoI_interpolated.attrs[attr] = input_global_attributes[attr]
    date_start = pd.to_datetime(VoI_interpolated.TIME.values.min()).strftime(TIMESTAMP_FORMAT)
    date_end = pd.to_datetime(VoI_interpolated.TIME.values.max()).strftime(TIMESTAMP_FORMAT)
    date_created = current_utc_timestamp()
    VoI_interpolated.attrs.update(global_attribute_dictionary)
    VoI_interpolated.attrs.update({
        'source_file':           input_file,
        'time_coverage_start':   date_start,
        'time_coverage_end':     date_end,
        'geospatial_vertical_min': min(depth_bins),
        'geospatial_vertical_max': max(depth_bins),
        'keywords':              ', '.join([VoI, 'DEPTH'] + ['HOURLY', 'GRIDDED']),
        'abstract':              global_attribute_dictionary['abstract'].format(VoI=VoI, site_code=site_code),
        'date_created': date_created,
        'history': input_global_attributes['history'] +' {date_created}: Gridded file created.'.format(
            date_created=date_created),
        'generating_code_version': __version__,
        'title':                 global_attribute_dictionary['title'].format(VoI=VoI,
                                                                    site_code=site_code,
                                                                    time_min=date_start,
                                                                    time_max=date_end,
                                                                    depth_min=min(depth_bins),
                                                                    depth_max = max(depth_bins))
    })
    github_comment = ('\nThis file was created using https://github.com/aodn/python-aodntools/blob/'
                      '{v}/aodntools/timeseries_products/{f}'.format(v=__version__, f=os.path.basename(__file__))
                      )
    VoI_interpolated.attrs['lineage'] += github_comment
    if download_url_prefix:
        VoI_interpolated.attrs['source_file_download'] = os.path.join(download_url_prefix, input_file)
    if opendap_url_prefix:
        VoI_interpolated.attrs['source_file_opendap'] = os.path.join(opendap_url_prefix, input_file)
    VoI_interpolated.attrs = sorted(VoI_interpolated.attrs.items())

    ## create the output file name and write the aggregated product as netCDF
    facility_code = TStools.get_facility_code(input_file)
    data_code = TStools.get_data_code(VoI) + 'Z'
    product_type = 'gridded-timeseries'
    file_version = 2
    ncout_filename = generate_netcdf_output_filename(nc=VoI_interpolated, facility_code=facility_code,
                                                             data_code=data_code, VoI=VoI,
                                                             site_code=site_code, product_type=product_type,
                                                             file_version=file_version)
    ncout_path = os.path.join(output_dir, ncout_filename)

    encoding = {'TIME': {'_FillValue': None,
                         'units': time_units,
                         'calendar': time_calendar,
                         'zlib': True,
                         'complevel': 5},
                VoI:    {'zlib': True,
                         'complevel': 5,
                         'dtype': np.dtype('float32')},
                VoI+'_count': {'dtype': np.dtype('int16'),
                               'zlib': True,
                               'complevel': 5},
                'DEPTH': {'dtype': np.dtype('float32'),
                          'zlib': True,
                          'complevel': 5},
                'LONGITUDE': {'_FillValue': False},
                'LATITUDE': {'_FillValue': False}}

    write_netCDF_aggfile(VoI_interpolated, ncout_path, encoding)

    return ncout_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gridded time series: interpolate ONE variable from ALL instruments from ALL deployments from ONE site into 1hr timestamps and fixed depth bins")
    parser.add_argument('-var', dest='var', help='name of the variable to concatenate. Like TEMP, PSAL', default='TEMP', required=False)
    parser.add_argument('-file', dest='filename', help='name of the Hourly Time Series Product file that contains the data', default=None, required=False)
    parser.add_argument('-depth_bins', dest='depth_bins', help='list of depth where the VoI will be interpolated', default=None, nargs='+', required=False)
    parser.add_argument('-max_separation', dest='max_separation', help='maximum difference between instruments to allow interpolation', default=50, required=False)
    parser.add_argument('-depth_bins_increment', dest='depth_bins_increment', help='increment in meters for the automatic generated depth bins', default=10, required=False)
    parser.add_argument('-indir', dest='input_dir', help='base path of input file. Default .', default='.',
                        required=False)
    parser.add_argument('-outdir', dest='output_dir', help='path where the result file will be written. Default .',
                        default='.', required=False)
    parser.add_argument('-config', dest='config_file', help='JSON configuration file', default=None, required=False)
    args = parser.parse_args()

    if args.config_file:
        with open(args.config_file) as ff:
            arguments = json.load(ff)
        VoI = arguments['var']
        depth_bins = arguments['depth_bins']
        depth_bins_increment = int(arguments['depth_bins_increment'])
        max_separation = int(arguments['max_separation'])
        input_dir = arguments.get('input_dir', '.')
        output_dir = arguments.get('output_dir', '.')
    else:
        VoI = args.var
        depth_bins = args.depth_bins
        depth_bins_increment = int(args.depth_bins_increment)
        max_separation = int(args.max_separation)
        input_dir = args.input_dir
        output_dir = args.output_dir

    file_name = args.filename

    print(grid_variable(input_file=file_name, VoI=VoI, depth_bins=depth_bins,
                        max_separation=int(max_separation), depth_bins_increment=int(depth_bins_increment),
                        input_dir=input_dir, output_dir=output_dir))
