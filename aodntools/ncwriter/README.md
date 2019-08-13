# ncwriter

This package provides a basic templating tool for creating netCDF files based on templates.

Templates can be provided as individual Python dictionaries, or as a JSON file, e.g. `template.json`:
```json
{
    "_dimensions": {
        "TIME": null
    },
    "_variables": {
        "TIME": {
            "_dimensions": ["TIME"],
            "_datatype": "float64"
        }
    },
    "title": "test dataset"
}
```
Basic usage
```python
import numpy as np
from aodntools.ncwriter import DatasetTemplate

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
template.variables["TIME"].update({
    "standard_name": "time",
    "units": "days since 1950-01-01 00:00:00 +00:00"
    })
template.variables["TEMP"] = {
    "_dimensions": ["TIME"],
    "_datatype": "float32",
    "standard_name": "sea_water_temperature",
    "units": "degC",
    "valid_min": 0.0,
    "valid_max": 42.00
}

# add variable values
template.variables['TIME']['_data'] = np.arange(10)
template.variables['TEMP']['_data'] = np.array([12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0])

# create netCDF file (dimensions set to 'null' in the template are adjusted to match data array sizes, if possible)
template.to_netcdf('example.nc')

```

See the [examples](../../examples) folder for more complete examples.
