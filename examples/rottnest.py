"""Example Python script to create a netCDF file using ncwriter

This script reads sample data from a CSV file (originally from IMOS/ANMN/NRS/NRSROT/CTD_timeseries/
IMOS_ANMN-NRS_CSTZ_20180406T080001Z_NRSROT_FV01_NRSROT-1804-SBE37SM-RS232-24_END-20180817T044501Z_C-20180820T010303Z.nc)
and writes it to a netCDF file.
"""

import os

import numpy as np
import pandas as pd
from netCDF4 import date2num

from aodntools.ncwriter import ImosTemplate, TIMESTAMP_FORMAT

EXAMPLES_PATH = os.path.dirname(__file__)
TEMPLATE_JSON = os.path.join(EXAMPLES_PATH, 'rottnest.json')
DATA_CSV = os.path.join(EXAMPLES_PATH, 'rottnest.csv')


# read data from CSV
df = pd.read_csv(DATA_CSV, parse_dates=['TIME'])

# create template
template = ImosTemplate.from_json(TEMPLATE_JSON)

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
template.add_extent_attributes()

# add creation date
template.add_date_created_attribute()
template.global_attributes['history'] = "{}: File created".format(template.date_created.strftime(TIMESTAMP_FORMAT))

# generate file name
outfile = 'rottnest.nc'

# write file
template.to_netcdf(outfile)
