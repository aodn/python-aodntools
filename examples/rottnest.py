"""Example Python script to create a netCDF file using ncwriter

This script reads sample data from a CSV file (originally from IMOS/ANMN/NRS/NRSROT/CTD_timeseries/
IMOS_ANMN-NRS_CSTZ_20180406T080001Z_NRSROT_FV01_NRSROT-1804-SBE37SM-RS232-24_END-20180817T044501Z_C-20180820T010303Z.nc)
and writes it to a netCDF file.
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
from netCDF4 import date2num

import ncwriter
from ncwriter.template import DatasetTemplate

IMOS_GLOBAL_JSON = os.path.join(ncwriter.__path__[0], 'imos_global.json')
EXAMPLES_PATH = os.path.dirname(__file__)
TEMPLATE_JSON = os.path.join(EXAMPLES_PATH, 'rottnest.json')
DATA_CSV = os.path.join(EXAMPLES_PATH, 'rottnest.csv')

TIMESTAMP_FORMAT = '%Y-%m-%dT%H:%M:%SZ'

# read data from CSV
df = pd.read_csv(DATA_CSV, parse_dates=['TIME'])

# create template
template = (DatasetTemplate.from_json(IMOS_GLOBAL_JSON) +
            DatasetTemplate.from_json(TEMPLATE_JSON)
            )

# update attributes
for att in ('site_code', 'platform_code', 'deployment_code', 'instrument_nominal_depth'):
    template.global_attributes[att] = df[att].unique()[0]

# add data
t_data = df['TIME'].dt.to_pydatetime()
t_template = template.variables['TIME']
t_template['_data'] = date2num(t_data, t_template['units'], t_template['calendar'])

template.variables['LATITUDE']['_data'] = df['LATITUDE'].unique()[0]
template.variables['LONGITUDE']['_data'] = df['LONGITUDE'].unique()[0]

for name, var in template.variables.items():
    if '_data' not in var:
        var['_data'] = df[name].values

# convert valid_min/max attributes to match variable type
# TODO: make this a template method
for name, var in template.variables.items():
    var_type = var['_datatype']
    for attr in ('valid_min', 'valid_max'):
        if attr in var:
            var[attr] = np.cast[var_type](var[attr])

# update range attributes
# TODO: make this a template method
template.global_attributes['time_coverage_start'] = t_data.min().strftime(TIMESTAMP_FORMAT)
template.global_attributes['time_coverage_end'] = t_data.max().strftime(TIMESTAMP_FORMAT)
for varname, shortname in [('LATITUDE', 'lat'), ('LONGITUDE', 'lon')]:
    for stat in (min, max):
        attname = "geospatial_{shortname}_{stat.__name__}".format(shortname=shortname, stat=stat)
        template.global_attributes[attname] = stat(df[varname])

# add creation date
# TODO: make this a template method
date_created = datetime.utcnow().strftime(TIMESTAMP_FORMAT)
template.global_attributes['date_created'] = date_created
template.global_attributes['history'] = "{}: File created".format(date_created)

# generate file name
outfile = 'rottnest.nc'

# write file
template.to_netcdf(outfile)
