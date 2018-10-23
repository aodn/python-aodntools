# aodn-netcdf-tools
[![build](https://travis-ci.org/aodn/aodn-netcdf-tools.png?branch=master)](https://travis-ci.org/aodn/aodn-netcdf-tools)

Repository for templates and code relating to generating standard NetCDF files for the Australia Ocean Data Network

## ncwrwiter

This module provides a basic templating tool for creating netCDF files based on templates.

Templates can be provided as individual Python dictionaries, or as a JSON file, e.g. `template.json`:
```json
{
    "dimensions": {
        "TIME": 0
    },
    "variables": {
        "TIME": {
            "dimensions": ["TIME"],
            "type": "float64"
        }
    },
    "global_attributes": {
        "title": "test dataset"
    }
}
```
Basic usage
```python
from ncwriter import DatasetTemplate

# create template from individual dictionaries
template = DatasetTemplate(dimensions=dims_dict, variables=var_dict, global_attributes=gatt_dict)
# OR from a template file
template = DatasetTemplate.from_json('template.json')

# add/update attributes
template.global_attributes.update({
    'comment': 'this was added later',
    'date_created': '2018-09-20T00:00:00Z'
    })

# add/update variables
template.variables["TIME"]["atttributes"] = {
    "standard_name": "time",
    "units": "days since 1950-01-01 00:00:00 +00:00"
}
template.variables["TEMP"] = {
    "dimensions": ["TIME"],
    "type": "float32",
    "attributes": {
        "standard_name": "sea_water_temperature",
        "units": "degC",
        "valid_min": 0.0,
        "valid_max": 42.00
    }
}

# add variable values (automatically updates size of corresponding dimensions)
template.variables['TIME']['values'] = np.arange(10)
template.variables['TEMP']['values'] = np.array([12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0])

# create netCDF file
template.to_netcdf('example.nc')

```
