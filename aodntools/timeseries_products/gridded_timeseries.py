import numpy as np
import bisect
import argparse
import os.path
import json
from datetime import datetime

import xarray as xr
import pandas as pd

import aodntools.timeseries_products.aggregated_timeseries as TStools


TEMPLATE_JSON = resource_filename(__name__, 'gridded_timeseries_template.json')

def get_site_code(file_name):
    """
    get site_code from the file name
    :param file_name: name of the file
    :return: site_code
    """
    return file_name.split("_")[4]


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

def sort_attributes(attributes):
    """

    :param attributes: dictionary of global attributes
    :return: sorted dictionary
    """
    return dict(sorted(attributes.items()))


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

    site_code = get_site_code(file_name)

    with xr.open_dataset(file_name) as nc:
        ## get global attributes
        global_attributes = nc.attrs

        ## in case no depth bins provided, create depth bins to the nearest rounded depth
        if not depth_bins:
            depth_bins = make_depth_bins(nc=nc, increment=depth_bins_increment)

        ## create empty containers
        time_0 = pd.to_datetime('1950-01-01T00:00:00')
        time_min = nc.TIME.values.min()
        depth_bin_len = len(depth_bins)
        # VoI_binned = xr.DataArray(np.full((depth_bin_len, 1), np.nan), coords=[depth_bins, [nc.TIME.min().values]],
        #                           dims=['DEPTH', 'TIME'])
        # VoI_ndepths = xr.DataArray(np.full(1, 0), coords=[[nc.TIME.min().values]], dims=['TIME'])
        # VoI_longitude = xr.DataArray(np.full(1, 0), coords=[[nc.TIME.min().values]], dims=['TIME'])
        # VoI_latitude = xr.DataArray(np.full(1, 0), coords=[[nc.TIME.min().values]], dims=['TIME'])
        VoI_binned = xr.DataArray(np.full((depth_bin_len, 1), np.nan), coords=[depth_bins, [time_0]],
                                  dims=['DEPTH', 'TIME'])
        VoI_ndepths = xr.DataArray(np.full(1, 0), coords=[[time_0]], dims=['TIME'])
        VoI_longitude = xr.DataArray(np.full(1, 0), coords=[[time_0]], dims=['TIME'])
        VoI_latitude = xr.DataArray(np.full(1, 0), coords=[[time_0]], dims=['TIME'])

        ## group nc by individual timestamps
        VoI_grouped = nc.groupby('TIME')

        for timestamp, (name, group) in enumerate(VoI_grouped):
            time = [name]
            n_depths = int(len(group.TEMP))
            latitude = float(group.LATITUDE.mean())
            longitude = float(group.LONGITUDE.mean())

            if n_depths >= 2:
                VoI_values = list(group[VoI].values)
                depth = list(group.DEPTH.values)

                ## check for max separation
                depth_mask = get_depth_mask(depth_bins=depth_bins, depths=depth, max_separation=max_separation)
                ## do the interpolation
                VoI_gridded = np.interp(depth_bins, depth, VoI_values, left=np.nan, right=np.nan)
                ## set masked depth bins to zero
                VoI_gridded = VoI_gridded * depth_mask
                VoI_gridded[VoI_gridded == 0] = np.nan

            else:
                VoI_gridded = np.full((depth_bin_len, 1), np.nan)

            VoI_temp = xr.DataArray(VoI_gridded.reshape(depth_bin_len, 1), coords=[depth_bins, time],
                                    dims=['DEPTH', 'TIME'])
            VoI_ndepths_temp = xr.DataArray([n_depths], coords=[time], dims=['TIME'])
            VoI_longitude_temp = xr.DataArray([longitude], coords=[time], dims=['TIME'])
            VoI_latitude_temp = xr.DataArray([latitude], coords=[time], dims=['TIME'])

            VoI_binned = xr.concat([VoI_binned, VoI_temp], dim='TIME')
            VoI_ndepths = xr.concat([VoI_ndepths, VoI_ndepths_temp], dim='TIME')
            VoI_longitude = xr.concat([VoI_longitude, VoI_longitude_temp], dim='TIME')
            VoI_latitude = xr.concat([VoI_latitude, VoI_latitude_temp], dim='TIME')

    VoI_interpolated = xr.Dataset({VoI: VoI_binned.astype(float),
                                   VoI + '_count': VoI_ndepths.astype(int),
                                   'LONGITUDE': VoI_longitude,
                                   'LATITUDE': VoI_latitude})
    ## drop the very first record as it is dummy
    VoI_interpolated = VoI_interpolated.where(VoI_interpolated.TIME >= time_min, drop=True)

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
    VoI_interpolated.attrs.update(global_attributes)
    add_attribute = {}
    add_attribute.update({'source_file': file_name})

    VoI_interpolated.attrs.update(TStools.set_globalattr(VoI_interpolated, TEMPLATE_JSON, VoI, site_code, add_attribute))
    VoI_interpolated.attrs['file_version'] = global_attribute_dictionary['file_version']
    VoI_interpolated.attrs['abstract'] = global_attribute_dictionary['abstract'].format(VoI=VoI, site_code=site_code)
    VoI_interpolated.attrs['title'] = global_attribute_dictionary['title'].format(VoI=VoI,
                                                                             site_code=site_code,
                                                                             time_min=pd.to_datetime(VoI_interpolated.TIME.values.min()).strftime(timeformat),
                                                                             time_max=pd.to_datetime(VoI_interpolated.TIME.values.max()).strftime(timeformat),
                                                                             depth_min=min(depth_bins),
                                                                             depth_max = max(depth_bins))
    VoI_interpolated.attrs = sort_attributes(VoI_interpolated.attrs)

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

    encoding = {'TIME': {'_FillValue': False,
                         'units': time_units,
                         'calendar': time_calendar},
                VoI+'_count': {'_FillValue': False,
                               'dtype': np.dtype('int')},
                'LONGITUDE': {'_FillValue': False},
                'LATITUDE': {'_FillValue': False}}

    TStools.write_netCDF_aggfile(VoI_interpolated, ncout_path, encoding)

    return ncout_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Gridded time series: interpolate ONE variable from ALL instruments from ALL deployments from ONE site into 1hr timestamps and fixed depth bins")
    parser.add_argument('-var', dest='varname', help='name of the variable to concatenate. Like TEMP, PSAL', required=True)
    parser.add_argument('-file', dest='filename', help='name of the Hourly Time Series Product file that contains the data', required=True)
    parser.add_argument('-depth_bins', dest='depth_bins', help='list of depth where the VoI will be interpolated', default=None, nargs='+', required=False)
    parser.add_argument('-max_separation', dest='max_separation', help='maximum difference between instruments to allow interpolation', default=50, required=False)
    parser.add_argument('-depth_bins_increment', dest='depth_bins_increment', help='increment in meters for the automatic generated depth bins', default=10, required=False)
    parser.add_argument('-path', dest='output_path', help='path where the result file will be written. Default ./', default='./', required=False)
    args = parser.parse_args()

    print(grid_variable(file_name=args.filename, VoI=args.varname, depth_bins=args.depth_bins,
                        max_separation=int(args.max_separation), depth_bins_increment=int(args.depth_bins_increment),
                        base_path=args.output_path))
