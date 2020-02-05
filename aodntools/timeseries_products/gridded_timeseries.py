import numpy as np
import bisect
import argparse
import os.path
import json
from datetime import datetime

import xarray as xr
import pandas as pd

from pkg_resources import resource_filename

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




## MAIN FUNCTION
def grid_variable(file_name, VoI, depth_bins=None, max_separation=50, depth_bins_increment=10,
                  base_path='./'):
    """
    Grid VoI into depth_bins.
    :param nc: hourly aggregated dataset with VoI, DEPTH and TIME only
    :param VoI: variable of interest (TEMP or PSAL)
    :param depth_bins: list of depth where to interpolate. if null list is provided it will be calculated from the data
    :param max_separation: max separation allowed for instruments
    :param depth_bins_increment: in case no depth bins provided this is the increment for the calculated bins
    :param base_path: path where the result file will be written
    :return: interpolated dataset
    """

    with xr.open_dataset(file_name) as nc_full:
        nc = nc_full[[VoI, 'TIME', 'DEPTH']]
        ## get lat/lon
        longitude_mean = nc_full.LONGITUDE.mean()
        latitude_mean = nc_full.LATITUDE.mean()

    site_code = nc.site_code
    ## get global attributes
    global_attributes = nc.attrs



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
    VoI_temp = xr.DataArray(np.full((depth_bin_len, 1), np.nan), coords=[depth_bins, [time_0]],
                              dims=['DEPTH', 'TIME'])
    VoI_ndepths = xr.DataArray(np.full(1, 0), coords=[[time_0]], dims=['TIME'])

    ## group nc by individual timestamps
    VoI_grouped = nc.groupby('TIME')


    for timestamp, group in VoI_grouped:
        time = [timestamp]
        n_depths = int(len(group[VoI]))

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

    VoI_interpolated = xr.Dataset({VoI: VoI_temp.astype(float),
                                   VoI + '_count': VoI_ndepths.astype(int)})

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
    variable_attributes = TStools.set_variableattr(varlist, variable_attribute_dictionary, add_variable_attribute)
    time_units = variable_attributes['TIME'].pop('units')
    time_calendar = variable_attributes['TIME'].pop('calendar')
    for variable in varlist:
        VoI_interpolated[variable].attrs = variable_attributes[variable]

    ## set global attributes
    timeformat = '%Y-%m-%dT%H:%M:%SZ'
    date_start = pd.to_datetime(VoI_interpolated.TIME.values.min()).strftime(timeformat)
    date_end = pd.to_datetime(VoI_interpolated.TIME.values.max()).strftime(timeformat)
    VoI_interpolated.attrs.update(global_attributes)
    VoI_interpolated.attrs.update({
        'file_version':          global_attribute_dictionary['file_version'],
        'source_file':           file_name,
        'featureType':           global_attribute_dictionary['featureType'],
        'time_coverage_start':   date_start,
        'time_coverage_end':     date_end,
        'geospatial_vertical_min': min(depth_bins),
        'geospatial_vertical_max': max(depth_bins),
        'keywords':              ', '.join([VoI, 'DEPTH'] + ['HOURLY', 'GRIDDED']),
        'abstract':              global_attribute_dictionary['abstract'].format(VoI=VoI, site_code=site_code),
        'history':               VoI_interpolated.attrs['history'] + ' {today}: Gridded file created.'.format(today=datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%S')),
        'lineage':               global_attribute_dictionary['lineage'],
        'title':                 global_attribute_dictionary['title'].format(VoI=VoI,
                                                                    site_code=site_code,
                                                                    time_min=date_start,
                                                                    time_max=date_end,
                                                                    depth_min=min(depth_bins),
                                                                    depth_max = max(depth_bins))})
    VoI_interpolated.attrs = sorted(VoI_interpolated.attrs.items())

    ## create the output file name and write the aggregated product as netCDF
    facility_code = TStools.get_facility_code(file_name)
    data_code = TStools.get_data_code(VoI) + 'Z'
    product_type = 'gridded-timeseries'
    file_version = 2
    ncout_filename = TStools.generate_netcdf_output_filename(nc=VoI_interpolated, facility_code=facility_code,
                                                             data_code=data_code, VoI=VoI,
                                                             site_code=site_code, product_type=product_type,
                                                             file_version=file_version)
    ncout_path = os.path.join(base_path, ncout_filename)

    encoding = {'TIME': {'_FillValue': None,
                         'units': time_units,
                         'calendar': time_calendar,
                         'zlib': True,
                         'complevel': 5},
                VoI:    {'zlib': True,
                         'complevel': 5},
                VoI+'_count': {'_FillValue': None,
                               'dtype': np.dtype('int'),
                               'zlib': True,
                               'complevel': 5},
                'DEPTH': {'dtype': np.dtype('float64'),
                          'zlib': True,
                          'complevel': 5},
                'LONGITUDE': {'_FillValue': False},
                'LATITUDE': {'_FillValue': False}}

    TStools.write_netCDF_aggfile(VoI_interpolated, ncout_path, encoding)

    return ncout_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gridded time series: interpolate ONE variable from ALL instruments from ALL deployments from ONE site into 1hr timestamps and fixed depth bins")
    parser.add_argument('-var', dest='var', help='name of the variable to concatenate. Like TEMP, PSAL', default='TEMP', required=False)
    parser.add_argument('-file', dest='filename', help='name of the Hourly Time Series Product file that contains the data', default=None, required=False)
    parser.add_argument('-depth_bins', dest='depth_bins', help='list of depth where the VoI will be interpolated', default=None, nargs='+', required=False)
    parser.add_argument('-max_separation', dest='max_separation', help='maximum difference between instruments to allow interpolation', default=50, required=False)
    parser.add_argument('-depth_bins_increment', dest='depth_bins_increment', help='increment in meters for the automatic generated depth bins', default=10, required=False)
    parser.add_argument('-path', dest='output_path', help='path where the result file will be written. Default ./', default='./', required=False)
    parser.add_argument('-config', dest='config_file', help='JSON configuration file', default=None, required=False)
    args = parser.parse_args()

    if args.config_file:
        with open(args.config_file) as ff:
            arguments = json.load(ff)
        VoI = arguments['var']
        depth_bins = arguments['depth_bins']
        depth_bins_increment = int(arguments['depth_bins_increment'])
        max_separation = int(arguments['max_separation'])
        output_path = arguments['output_path']
    else:
        VoI = args.var
        depth_bins = args.depth_bins
        depth_bins_increment = int(args.depth_bins_increment)
        max_separation = int(args.max_separation)
        output_path = args.output_path

    file_name = args.filename

    print(grid_variable(file_name=file_name, VoI=VoI, depth_bins=depth_bins,
                        max_separation=int(max_separation), depth_bins_increment=int(depth_bins_increment),
                        base_path=output_path))
