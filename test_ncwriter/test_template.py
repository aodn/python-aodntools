import json
import os
import shutil
import tempfile
import unittest
from collections import OrderedDict

import numpy as np
from netCDF4 import Dataset

from ncwriter import DatasetTemplate, ValidationError

TEST_ROOT = os.path.dirname(__file__)
TEMPLATE_JSON = os.path.join(TEST_ROOT, 'template1.json')
TEMPLATE_PARTIAL_JSON = os.path.join(TEST_ROOT, 'template_partial.json')


class TestDatasetTemplate(unittest.TestCase):
    with open(TEMPLATE_JSON) as t:
        template_dict = json.load(t, object_pairs_hook=OrderedDict)
    dimensions = template_dict['dimensions']
    variables = template_dict['variables']
    global_attributes = template_dict['global_attributes']
    values1 = np.array([1], dtype=np.float32)
    values10 = np.arange(10, dtype=np.float32)

    @property
    def temp_dir(self):
        if not hasattr(self, '_temp_dir'):
            self._temp_dir = tempfile.mkdtemp(prefix=self.__class__.__name__)
        return self._temp_dir

    @property
    def temp_nc_file(self):
        if not hasattr(self, '_temp_nc_file'):
            with tempfile.NamedTemporaryFile(suffix='.nc', prefix=self.__class__.__name__, dir=self.temp_dir) as f:
                pass
            self._temp_nc_file = f.name
        return self._temp_nc_file

    def tearDown(self):
        if hasattr(self, '_temp_dir'):
            shutil.rmtree(self._temp_dir)

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

    def test_init_from_dicts_validation(self):
        with self.assertRaises(ValidationError):
            DatasetTemplate(dimensions='X')
        with self.assertRaises(ValidationError):
            DatasetTemplate(dimensions={'TIME': -1})

        with self.assertRaises(ValidationError):
            DatasetTemplate(variables='TEMP')
        with self.assertRaises(ValidationError):
            DatasetTemplate(variables={'_TEMP': {}})

        with self.assertRaises(ValidationError):
            DatasetTemplate(global_attributes='title')
        with self.assertRaises(ValidationError):
            DatasetTemplate(global_attributes={'title': None})

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
        template.variables['TEMP']['dimensions'] = ['TIME', 'DEPTH']
        self.assertEqual(['TIME', 'DEPTH'], template.variables['TEMP']['dimensions'])

    def test_add_variable_attributes(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.variables['TEMP']['attributes'].update([('units', 'Kelvin'),
                                                         ('comment', 'ok')
                                                         ])
        self.assertEqual(OrderedDict([('standard_name', 'sea_water_temperature'),
                                      ('units', 'Kelvin'),
                                      ('comment', 'ok')
                                      ]),
                         template.variables['TEMP']['attributes']
                         )

    def test_set_variable_values(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        template.variables['TEMP']['values'] = self.values10
        self.assertTrue(all(template.variables['TEMP']['values'] == self.values10))

    def test_create_empty_file(self):
        template = DatasetTemplate()
        template.to_netcdf(self.temp_nc_file)
        dataset = Dataset(self.temp_nc_file)

    def test_create_file(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        template.variables['TIME']['values'] = self.values10
        template.variables['DEPTH']['values'] = self.values1
        template.variables['TEMP']['values'] = self.values10.reshape((10, 1))
        template.to_netcdf(self.temp_nc_file)

        dataset = Dataset(self.temp_nc_file)

        expected_dimensions = OrderedDict([
            ('TIME', len(self.values10)),
            ('DEPTH', len(self.values1))
        ])
        ds_dimensions = OrderedDict((k, v.size) for k, v in dataset.dimensions.iteritems())
        self.assertEqual(expected_dimensions, ds_dimensions)

        for vname, vdict in self.variables.iteritems():
            ds_var = dataset[vname]
            self.assertEqual(vdict['dimensions'], list(ds_var.dimensions))
            self.assertEqual(vdict['type'], ds_var.dtype)
            ds_var_attr = OrderedDict((k, ds_var.getncattr(k)) for k in ds_var.ncattrs())
            if vdict.get('attributes') is None:
                self.assertEqual({}, ds_var_attr)
            else:
                self.assertEqual(vdict['attributes'], ds_var_attr)

        self.assertTrue(all(dataset['TIME'] == self.values10))
        self.assertTrue(all(dataset['DEPTH'] == self.values1))
        self.assertTrue(all(dataset['TEMP'] == self.values10.reshape(10, 1)))

        ds_global_attributes = OrderedDict((k, dataset.getncattr(k)) for k in dataset.ncattrs())
        self.assertEqual(self.global_attributes, ds_global_attributes)

    def test_fill_values(self):
        template = DatasetTemplate(dimensions={'X': 10})
        template.variables['X'] = {'dimensions': ['X'],
                                   'type': 'float32',
                                   'attributes': {'_FillValue': -999.}
                                   }
        x = np.array([-999., -999., -999., -999., -999., 1., 2., 3., 4., 5])
        template.variables['X']['values'] = x
        template.to_netcdf(self.temp_nc_file)

        dataset = Dataset(self.temp_nc_file)
        dsx = dataset.variables['X']
        self.assertEqual(-999., dsx._FillValue)
        self.assertIsInstance(dsx[:], np.ma.MaskedArray)
        self.assertTrue(dsx[:5].mask.all())
        self.assertTrue((dsx[5:] == x[5:]).all())

    def test_fill_values_from_masked_array(self):
        template = DatasetTemplate(dimensions={'X': 10})
        template.variables['X'] = {'dimensions': ['X'],
                                   'type': 'float32',
                                   'attributes': {'_FillValue': -999.}
                                   }
        x = np.array([-4, -3, -2, -1, 0, 1., 2., 3., 4., 5])
        template.variables['X']['values'] = np.ma.masked_array(x,
                                                               mask=[True, True, True, True, True,
                                                                     False, False, False, False, False]
                                                               )
        template.to_netcdf(self.temp_nc_file)

        dataset = Dataset(self.temp_nc_file)
        dsx = dataset.variables['X']
        self.assertEqual(-999., dsx._FillValue)
        self.assertIsInstance(dsx[:], np.ma.MaskedArray)
        self.assertTrue(dsx[:5].mask.all())
        self.assertTrue((dsx[5:] == x[5:]).all())

# TODO: add data from multiple numpy arrays
# e.g. template.add_data(TIME=time_values, TEMP=temp_values, PRES=pres_values)
# TODO: add data from Pandas dataframe (later...)
# e.g. template.add_data(dataframe)

# TODO: create netCDF file with auto-generated file name according to IMOS conventions
# e.g. template.create()


if __name__ == '__main__':
    unittest.main()
