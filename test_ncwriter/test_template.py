from collections import OrderedDict
import numpy as np
import unittest

from ncwriter import DatasetTemplate


class TestDatasetTemplate(unittest.TestCase):
    dimensions = OrderedDict([
        ('time', 0),
        ('depth', 0)
    ])
    time = OrderedDict([
        ('dims', ['time']),
        ('type', np.float64),
        ('attr', None)
    ])
    temp = OrderedDict([
        ('dims', ['time', 'depth']),
        ('type', np.float32),
        ('attr', OrderedDict([('standard_name', 'sea_water_temperature'), ('units', 'degC')]))
    ])
    variables = OrderedDict([
        ('TIME', time),
        ('TEMP', temp)
    ])
    global_attributes = OrderedDict([
        ('title', 'test dataset'),
        ('abstract', 'This is a dataset used for testing'),
        ('Conventions', 'CF-1.6,IMOS-1.4')
    ])

    def test_init_empty(self):
        template = DatasetTemplate()
        self.assertEqual({}, template.dimensions)
        self.assertEqual({}, template.variables)
        self.assertEqual({}, template.global_attributes)

    def test_init_from_dicts(self):
        template = DatasetTemplate(dimensions=self.dimensions,
                                   variables=self.variables,
                                   global_attributes=self.global_attributes)
        self.assertEqual(self.dimensions, template.dimensions)
        self.assertEqual(self.variables, template.variables)
        self.assertEqual(self.global_attributes, template.global_attributes)

# TODO: create template from json file
# e.g. DatasetTemplate.from_json(path='template.json')
# TODO: create template from other formats (later...)

# TODO: add global attributes
# e.g. template.title = 'Test dataset'

# TODO: add dimensions
# e.g. template.dimensions['TIME'] = 100

# TODO: add variables
# TODO: add variable attributes
# e.g. template.variables['PRES']['units'] = 'dbar'

# TODO: add data from numpy arrays
# e.g. template.add_data(TIME=time_values, TEMP=temp_values, PRES=pres_values)
# TODO: add data from Pandas dataframe (later...)
# e.g. template.add_data(dataframe)

# TODO: create netCDF file
# e.g. template.create(filename)   # user-specified file name
# TODO: create netCDF file with auto-generated file name according to IMOS conventions
# e.g. template.create()


if __name__ == '__main__':
    unittest.main()
