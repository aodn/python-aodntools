#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geoserverCatalog.py
Collect files names from the AODN geoserver according to several conditions
Output a list of urls
"""

from __future__ import print_function

import sys
import argparse
from datetime import datetime

import pandas as pd


def args():
    parser = argparse.ArgumentParser(description="Get a list of urls from the AODN geoserver")
    parser.add_argument('-var', dest='varname', help='name of the variable of interest, like TEMP', default=None, required=False)
    parser.add_argument('-site', dest='site', help='site code, like NRMMAI',  type=str, default=None, required=False)
    parser.add_argument('-ft', dest='featuretype', help='feature type, like timeseries', default=None, required=False)
    parser.add_argument('-fv', dest='fileversion', help='file version, like 1', default=None, type=int, required=False)
    parser.add_argument('-ts', dest='timestart', help='start time like 2015-12-01', default=None, type=str, required=False)
    parser.add_argument('-te', dest='timeend', help='end time like 2018-06-30', type=str, default=None, required=False)
    parser.add_argument('-dc', dest='datacategory', help='data category like Temperature', type=str, default=None, required=False)
    parser.add_argument('-realtime', dest='realtime', help='yes or no. If absent, all modes will be retrieved', type=str, default=None, required=False)
    parser.add_argument('-rm', dest='filterout', help='regex to filter out the url list. Case sensitive', type=str, default=None, required=False)

    vargs = parser.parse_args()
    return(vargs)


def get_moorings_urls(varname=None, site=None, featuretype=None, fileversion=None, datacategory=None, realtime=None, timestart=None, timeend=None, filterout=None):
    """
    get the urls from the geoserver moorings_all_map collection
    based on user defined filters
    """

    WEBROOT = 'http://thredds.aodn.org.au/thredds/dodsC/'

    if realtime:
        if realtime.lower() == "yes":
            url = "http://geoserver-123.aodn.org.au/geoserver/ows?typeName=moorings_all_map&SERVICE=WFS&REQUEST=GetFeature&VERSION=1.0.0&outputFormat=csv&CQL_FILTER=(realtime=TRUE)"
        elif realtime.lower() == "no":
            url = "http://geoserver-123.aodn.org.au/geoserver/ows?typeName=moorings_all_map&SERVICE=WFS&REQUEST=GetFeature&VERSION=1.0.0&outputFormat=csv&CQL_FILTER=(realtime=FALSE)"
        else:
            raise ValueError('ERROR: realtime %s is not valid' % realtime)
    else:
        url = "http://geoserver-123.aodn.org.au/geoserver/ows?typeName=moorings_all_map&SERVICE=WFS&REQUEST=GetFeature&VERSION=1.0.0&outputFormat=csv"

    df = pd.read_csv(url)
    df = df.sort_values(by='time_coverage_start')
    criteria_all = df.url != None

    if varname:
        separator = ', '
        varnames_all = set(separator.join(list(df.variables)).split(', '))
        if varname in varnames_all:
            criteria_all = criteria_all & df.variables.str.contains('.*\\b'+varname+'\\b.*', regex=True)
        else:
            raise ValueError('ERROR: %s not a valid variable name' % varname)

    if site:
        site_all = list(df.site_code.unique())
        if site in site_all:
            criteria_all = criteria_all & df.site_code.str.contains(site, regex=False)
        else:
            raise ValueError('ERROR: %s is not a valid site code' % site)

    if featuretype:
        if featuretype in ["timeseries", "profile", "timeseriesprofile"]:
            criteria_all = criteria_all & (df.feature_type.str.lower() == featuretype.lower())
        else:
            raise ValueError('ERROR: %s is not a valid feature type' % featuretype)

    if datacategory:
        datacategory_all = list(df.data_category.str.lower().unique())
        if datacategory.lower() in datacategory_all:
            criteria_all = criteria_all & (df.data_category.str.lower() == datacategory.lower())
        else:
            raise ValueError('ERROR: %s is not a valid data category' % datacategory)

    if fileversion is not None:
        if fileversion in [0, 1, 2]:
            criteria_all = criteria_all & (df.file_version == fileversion)
        else:
            raise ValueError('ERROR: %s is not a valid file version' % fileversion)

    if timestart:
        try:
            criteria_all = criteria_all & (pd.to_datetime(df.time_coverage_end) >= datetime.strptime(timestart, '%Y-%m-%d'))
        except ValueError:
            raise ValueError('ERROR: invalid start date.')

    if timeend:
        try:
            criteria_all = criteria_all & (pd.to_datetime(df.time_coverage_start) <=  datetime.strptime(timeend, '%Y-%m-%d'))
        except ValueError:
            raise ValueError('ERROR: invalid end date.')

    if filterout is not None:
        criteria_all = criteria_all & (~df.url.str.contains(filterout, regex=True))


    return list(WEBROOT + df.url[criteria_all])


if __name__ == "__main__":
    vargs = args()
    urls = get_moorings_urls(**vars(vargs))
    print(*urls, sep='\n')
