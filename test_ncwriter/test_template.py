from collections import OrderedDict
import json
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
