from collections import OrderedDict
import json
import numpy as np
import os
import unittest

from ncwriter import DatasetTemplate


TEST_ROOT = os.path.dirname(__file__)
TEMPLATE_JSON = os.path.join(TEST_ROOT, 'template1.json')
TEMPLATE_PARTIAL_JSON = os.path.join(TEST_ROOT, 'template_partial.json')


class TestDatasetTemplate(unittest.TestCase):
    with open(TEMPLATE_JSON) as t:
        template_dict = json.load(t, object_pairs_hook=OrderedDict)
    dimensions = template_dict['dimensions']
    variables = template_dict['variables']
    global_attributes = template_dict['global_attributes']

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

    def test_init_from_json(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        self.assertEqual(self.dimensions, template.dimensions)
        self.assertEqual(self.variables, template.variables)
        self.assertEqual(self.global_attributes, template.global_attributes)

    def test_init_from_partial_template(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        with open(TEMPLATE_PARTIAL_JSON) as t:
            tdict = json.load(t, object_pairs_hook=OrderedDict)
        self.assertEqual({}, template.dimensions)
        self.assertEqual(tdict['variables'], template.variables)
        self.assertEqual(tdict['global_attributes'], template.global_attributes)

    # TODO: def test_json_validation(self):

    # TODO: create template from other formats (later...)

    def test_add_global_attributes(self):
        template = DatasetTemplate()
        template.global_attributes.update(self.global_attributes)
        self.assertEqual(self.global_attributes, template.global_attributes)

    def test_add_dimensions(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.dimensions['TIME'] = 100
        template.dimensions['DEPTH'] = 10
        self.assertEqual(OrderedDict([('TIME', 100), ('DEPTH', 10)]), template.dimensions)

    def test_update_dimensions(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        template.dimensions['TIME'] = 100
        template.dimensions['DEPTH'] = 10
        self.assertDictContainsSubset(OrderedDict([('TIME', 100), ('DEPTH', 10)]), template.dimensions)

    def test_add_variables(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.variables['TIME'] = self.variables['TIME']
        self.assertEqual(['TEMP', 'TIME'], template.variables.keys())
        self.assertEqual(self.variables['TIME'], template.variables['TIME'])

    def test_add_variable_dimensions(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.variables['TEMP']['dims'] = ['TIME', 'DEPTH']
        self.assertEqual(['TIME', 'DEPTH'], template.variables['TEMP']['dims'])

    def test_add_variable_attributes(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.variables['TEMP']['attr'].update([('units', 'Kelvin'),
                                                   ('comment', 'ok')
                                                   ])
        self.assertEqual(OrderedDict([('standard_name', 'sea_water_temperature'),
                                      ('units', 'Kelvin'),
                                      ('comment', 'ok')
                                      ]),
                         template.variables['TEMP']['attr']
                         )

    def test_set_variable_values(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        temp_val = np.arange(10, dtype=np.float32)
        template.variables['TEMP']['values'] = temp_val
        self.assertTrue(all(template.variables['TEMP']['values'] == temp_val))

# TODO: add data from multiple numpy arrays
# e.g. template.add_data(TIME=time_values, TEMP=temp_values, PRES=pres_values)
# TODO: add data from Pandas dataframe (later...)
# e.g. template.add_data(dataframe)

# TODO: create netCDF file
# e.g. template.create(filename)   # user-specified file name
# TODO: create netCDF file with auto-generated file name according to IMOS conventions
# e.g. template.create()


if __name__ == '__main__':
    unittest.main()
