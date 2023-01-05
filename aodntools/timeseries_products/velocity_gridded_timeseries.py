import numpy as np
import bisect
import argparse
import os.path
import json
from datetime import datetime

import xarray as xr
import pandas as pd

from pkg_resources import resource_filename

from aodntools import __version__
import aodntools.timeseries_products.aggregated_timeseries as TStools

import time as tt

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

    file_timeformat = '%Y%m%d'

    if '_' in VoI:
        VoI = VoI.replace('_', '-')
    t_start = pd.to_datetime(nc.TIME.min().values).strftime(file_timeformat)
    t_end = pd.to_datetime(nc.TIME.max().values).strftime(file_timeformat)

    output_name = '_'.join(['IMOS', facility_code, data_code, t_start, site_code, ('FV0'+str(file_version)), (VoI+"-"+product_type), ('END-'+ t_end), 'C-' + datetime.utcnow().strftime(file_timeformat)]) + '.nc'

    return output_name



## MAIN FUNCTION
def grid_variable(input_file, site_code, depth_bins=None, max_separation=50, depth_bins_increment=10,
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
        nc = nc_full[['VCUR','UCUR','TIME', 'DEPTH']]
        ## get lat/lon
        longitude_mean = nc_full.LONGITUDE.mean()
        latitude_mean = nc_full.LATITUDE.mean()
    ## get global attributes
    input_global_attributes = nc.attrs
    node = input_file[33:36];

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
    V_temp = xr.DataArray(np.full((depth_bin_len, 1), np.nan, dtype=np.float32), coords=[depth_bins, [time_0]],
                              dims=['DEPTH', 'TIME'])
    U_temp = xr.DataArray(np.full((depth_bin_len, 1), np.nan, dtype=np.float32), coords=[depth_bins, [time_0]],
                              dims=['DEPTH', 'TIME'])
    VoI_ndepths = xr.DataArray(np.full(1, 0, dtype='int'), coords=[[time_0]], dims=['TIME'])

    ## group nc by individual timestamps
    VoI_grouped = nc.groupby('TIME')

    t = tt.time()
    for timestamp, group in VoI_grouped:
        time = [timestamp]
        n_depths = np.array(len(group['VCUR']), dtype='int')
        
        if n_depths >= 2:
            V_values = group['VCUR'].values
            U_values = group['UCUR'].values
            D_values = group.DEPTH.values
            # remove nans if present
            if np.sum(np.logical_and(np.isnan(D_values), np.isnan(V_values))) > 0:
                c = np.logical_and(np.isfinite(D_values), np.isfinite(V_values)).transpose()
                V_values = list(np.array(V_values)[c])
                U_values = list(np.array(U_values)[c])
                D_values = list(np.array(D_values)[c]) 
                
            if len(V_values) != 0:
                ## sort depths
                depth, V_values = sort_depths(D_values, V_values)
                _, U_values = sort_depths(D_values, U_values)
                ## check for max separation
                depth_mask = get_depth_mask(depth_bins=depth_bins, depths=depth, max_separation=max_separation)
                ## do the interpolation
                interpolated_V = np.interp(depth_bins, depth, V_values, left=np.nan, right=np.nan)
                interpolated_U = np.interp(depth_bins, depth, U_values, left=np.nan, right=np.nan)
                ## set masked depth bins to zero
                interpolated_V = interpolated_V * depth_mask
                interpolated_V[interpolated_V == 0] = np.nan
                interpolated_U = interpolated_U * depth_mask
                interpolated_U[interpolated_U == 0] = np.nan
            else:
                interpolated_V = np.ones(np.size(depth_bins))*np.nan
                interpolated_U = np.ones(np.size(depth_bins))*np.nan
        else:
            interpolated_V = np.full((depth_bin_len, 1), np.nan)
            interpolated_U = np.full((depth_bin_len, 1), np.nan)

        V_temp_tmp = xr.DataArray(interpolated_V.reshape(depth_bin_len, 1), coords=[depth_bins, time],
                                dims=['DEPTH', 'TIME'])
        U_temp_tmp = xr.DataArray(interpolated_U.reshape(depth_bin_len, 1), coords=[depth_bins, time],
                                dims=['DEPTH', 'TIME'])
        VoI_ndepths_tmp = xr.DataArray([n_depths], coords=[time], dims=['TIME'])

        ## concatenate the interpolated values
        V_temp = xr.concat([V_temp, V_temp_tmp], dim='TIME')
        U_temp = xr.concat([U_temp, U_temp_tmp], dim='TIME')
        VoI_ndepths = xr.concat([VoI_ndepths, VoI_ndepths_tmp], dim='TIME')
        
    print(str(np.round(((tt.time() - t)/60),2)) + ' mins elapsed')

    V_interpolated = xr.Dataset({'VCUR': V_temp.astype(np.float32),
                                   'VCUR' + '_count': VoI_ndepths.astype('int16')})
    U_interpolated = xr.Dataset({'UCUR': U_temp.astype(np.float32),
                                   'UCUR' + '_count': VoI_ndepths.astype('int16')})

    ## drop the very first record as it is dummy
    V_interpolated = V_interpolated.where(V_interpolated.TIME >= time_min, drop=True)
    U_interpolated = U_interpolated.where(U_interpolated.TIME >= time_min, drop=True)

    ## Add lat/lon as scalar variables
    V_interpolated = V_interpolated.assign(LONGITUDE = longitude_mean,
                                               LATITUDE = latitude_mean)
    U_interpolated = U_interpolated.assign(LONGITUDE = longitude_mean,
                                               LATITUDE = latitude_mean)

    ## transpose dimensions to make CF compliant
    V_interpolated = V_interpolated.transpose('TIME', 'DEPTH')
    U_interpolated = U_interpolated.transpose('TIME', 'DEPTH')
    
    # concatenate V and U data sets
    griddedUV = xr.merge([V_interpolated,U_interpolated])
    
    # add global attributes 
    griddedUV.attrs = input_global_attributes
    # adapt title, lineage and abstract
    griddedUV.attrs['title'] = ('Gridded Time Series Product: VCUR, UCUR, WCUR interpolated at ' + site_code +
                                ' to fixed target depths at 1-hour time intervals, between ' + 
                                griddedUV.time_coverage_start + ' and ' + griddedUV.time_coverage_end + 
                                ' and ' + str(np.round(griddedUV.geospatial_vertical_min,1)) + ' and ' + 
                                str(np.round(griddedUV.geospatial_vertical_max,1)) + ' meters.')
    griddedUV.attrs['abstract'] = ('Gridded Time Series Product: This file contains VCUR, UCUR, WCUR readings' + 
                                   ' from all instruments deployed at the ' + site_code +' mooring site. The' +
                                   ' source of the values is the Hourly Time Series Product where TIME is ' + 
                                   'fixed to 1-hour interval. The variable values are interpolated to a' +
                                   ' fixed target depths using a linear interpolation between consecutive' + 
                                   ' existing depths. Only values flagged as 1 or 2 are used in the' + 
                                   ' interpolation.')
    griddedUV.attrs['lineage'] = ('The Variables of Interest VCUR, UCUR, WCUR are produced by sequentially ' + 
                                  'interpolating linearly the individual values at pre-defined target depths ' +
                                  'using the 1-hour binned data from the Hourly Time Series Product. The number ' +
                                  'of instruments used for the interpolation at each timestamp is recorded ' +
                                  'in the VoI_count variable. If less that two instrument readings are ' +
                                  'present in the timestamp, fill values are stored in all depth bins. The ' + 
                                  'resulting variable has dimensions TIME and DEPTH.')
    # add variable attributes
    griddedUV.TIME.attrs = nc.TIME.attrs
    griddedUV.DEPTH.attrs = nc.DEPTH.attrs
    griddedUV.VCUR.attrs = nc.VCUR.attrs
    griddedUV.UCUR.attrs = nc.UCUR.attrs
    
    # save file
    filename = ('IMOS_ANMN-' +  node + '_Z_' + nc.time_coverage_start[0:4] + \
                nc.time_coverage_start[5:7] + nc.time_coverage_start[8:10] + '_' + site_code + \
                '_FV01_velocity-gridded-timeseries_END-' + nc.time_coverage_end[0:4] + \
                nc.time_coverage_end[5:7] + nc.time_coverage_end[8:10] + '_C-' + \
                datetime.now().strftime("%Y%m%d") + '.nc')
    griddedUV.to_netcdf(output_dir + filename)

   