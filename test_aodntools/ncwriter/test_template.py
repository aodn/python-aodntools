#!/usr/bin/env python3

import json
import os
import re
import shutil
import tempfile
import unittest
from collections import OrderedDict

import numpy as np
from netCDF4 import Dataset

from test_aodntools.base_test import BaseTestCase
from aodntools.ncwriter import DatasetTemplate, ValidationError, metadata_attributes, special_attributes

TEST_ROOT = os.path.dirname(__file__)
TEMPLATE_JSON = os.path.join(TEST_ROOT, 'template1.json')
TEMPLATE_PARTIAL_JSON = os.path.join(TEST_ROOT, 'template_partial.json')
BAD_JSON = os.path.join(TEST_ROOT, 'bad.json')


class TestUtils(BaseTestCase):
    def test_metadata_attributes(self):
        self.assertEqual({}, metadata_attributes({}))
        self.assertEqual({}, metadata_attributes({'_dimensions': {}, '_fill_value': -999}))
        self.assertEqual({'title': 'Title'},
                         metadata_attributes({'title': 'Title'})
                         )
        self.assertEqual({'title': 'Title'},
                         metadata_attributes({'title': 'Title', '_fill_value': -999})
                         )
        self.assertIsInstance(metadata_attributes(OrderedDict()), OrderedDict)

    def test_special_attributes(self):
        self.assertEqual({}, special_attributes({}))
        self.assertEqual({}, special_attributes({'title': 'Title'}))
        self.assertEqual({'dimensions': {}, 'fill_value': -999},
                         special_attributes({'_dimensions': {}, '_fill_value': -999}))
        self.assertEqual({'fill_value': -999},
                         special_attributes({'title': 'Title', '_fill_value': -999}))


class TemplateTestCase(unittest.TestCase):
    with open(TEMPLATE_JSON) as t:
        template_dict = json.load(t, object_pairs_hook=OrderedDict)
    dimensions = template_dict['_dimensions']
    variables = template_dict['_variables']
    global_attributes = metadata_attributes(template_dict)
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


class TestDatasetTemplate(TemplateTestCase):
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

    def test_invalid_json(self):
        error_pattern = r"invalid JSON file '{}'".format(re.escape(BAD_JSON))
        self.assertRaisesRegexp(ValueError, error_pattern, DatasetTemplate.from_json, BAD_JSON)

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
        self.assertEqual(tdict['_variables'], template.variables)
        self.assertEqual(metadata_attributes(tdict), template.global_attributes)

    def test_add_method(self):
        template1 = DatasetTemplate(dimensions={'ONE': 1},
                                    variables={'X': {'_dimensions': ['ONE'], '_datatype': 'float32'},
                                               'Y': {'_dimensions': ['ONE'], '_datatype': 'float32'}
                                               },
                                    global_attributes={'title': 'First template', 'comment': 'one'}
                                    )
        template2 = DatasetTemplate(dimensions={'TWO': 2},
                                    variables={'Y': {'_dimensions': ['TWO'], 'comment': 'updated'},
                                               'Z': {'name': 'new'}
                                               },
                                    global_attributes={'title': 'Second template', 'version': 2}
                                    )
        template = template1 + template2

        self.assertEqual({'ONE': 1, 'TWO': 2}, template.dimensions)
        self.assertEqual({'title': 'Second template', 'comment': 'one', 'version': 2}, template.global_attributes)

        self.assertSetEqual({'X', 'Y', 'Z'}, set(template.variables.keys()))
        self.assertEqual({'_dimensions': ['ONE'], '_datatype': 'float32'}, template.variables['X'])
        self.assertEqual({'_dimensions': ['TWO'], '_datatype': 'float32', 'comment': 'updated'},
                         template.variables['Y'])
        self.assertEqual({'name': 'new'}, template.variables['Z'])

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

    def test_change_dimensions(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        template.dimensions['TIME'] = 100
        template.dimensions['DEPTH'] = 10
        self.assertEqual(OrderedDict([('TIME', 100), ('DEPTH', 10)]), template.dimensions)

    def test_add_variables(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.variables['TIME'] = self.variables['TIME']
        self.assertEqual({'TEMP', 'TIME'}, set(template.variables.keys()))
        self.assertEqual(self.variables['TIME'], template.variables['TIME'])

    def test_add_variable_dimensions(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.variables['TEMP']['_dimensions'] = ['TIME', 'DEPTH']
        self.assertEqual(['TIME', 'DEPTH'], template.variables['TEMP']['_dimensions'])

    def test_add_variable_attributes(self):
        template = DatasetTemplate.from_json(TEMPLATE_PARTIAL_JSON)
        template.variables['TEMP'].update([('units', 'Kelvin'),
                                           ('comment', 'ok')
                                           ])
        self.assertEqual(OrderedDict([('standard_name', 'sea_water_temperature'),
                                      ('units', 'Kelvin'),
                                      ('comment', 'ok')
                                      ]),
                         template.variables['TEMP']
                         )

    def test_set_variable_values(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        template.variables['TEMP']['_data'] = self.values10
        self.assertTrue(all(template.variables['TEMP']['_data'] == self.values10))

    def test_create_empty_file(self):
        template = DatasetTemplate()
        template.to_netcdf(self.temp_nc_file)
        dataset = Dataset(self.temp_nc_file)

    def test_create_empty_variable(self):
        template = DatasetTemplate(dimensions={'X': 10})
        template.variables['X'] = {'_dimensions': ['X'], '_datatype': 'float32'}
        self.assertRaises(ValidationError, template.to_netcdf, self.temp_nc_file)  # not providing '_data' is an error

        del self._temp_nc_file  # Get a new temp file
        template.variables['X']['_data'] = None  # This is ok, it's a shortcut for all fill values
        template.to_netcdf(self.temp_nc_file)

        dataset = Dataset(self.temp_nc_file)
        dataset.set_auto_mask(True)
        dsx = dataset.variables['X']
        self.assertIsInstance(dsx[:], np.ma.MaskedArray)
        self.assertTrue(dsx[:].mask.all())

    def test_create_file(self):
        template = DatasetTemplate.from_json(TEMPLATE_JSON)
        template.variables['TIME']['_data'] = self.values10
        template.variables['DEPTH']['_data'] = self.values1
        template.variables['TEMP']['_data'] = self.values10.reshape((10, 1))
        template.to_netcdf(self.temp_nc_file)

        dataset = Dataset(self.temp_nc_file)

        expected_dimensions = OrderedDict([
            ('TIME', len(self.values10)),
            ('DEPTH', len(self.values1))
        ])
        ds_dimensions = OrderedDict((k, v.size) for k, v in dataset.dimensions.items())
        self.assertEqual(expected_dimensions, ds_dimensions)

        for vname, vdict in self.variables.items():
            ds_var = dataset[vname]
            self.assertEqual(vdict['_dimensions'], list(ds_var.dimensions))
            self.assertEqual(vdict['_datatype'], ds_var.dtype)
            ds_var_attr = OrderedDict((k, ds_var.getncattr(k)) for k in ds_var.ncattrs())
            self.assertEqual(metadata_attributes(vdict), ds_var_attr)

        self.assertTrue(all(dataset['TIME'] == self.values10))
        self.assertTrue(all(dataset['DEPTH'] == self.values1))
        self.assertTrue(all(dataset['TEMP'] == self.values10.reshape(10, 1)))

        ds_global_attributes = OrderedDict((k, dataset.getncattr(k)) for k in dataset.ncattrs())
        self.assertEqual(self.global_attributes, ds_global_attributes)

    def test_close_file_on_exception(self):
        template = DatasetTemplate(variables={'Z': {}})
        self.assertIsNone(template.ncobj)
        self.assertRaises(ValidationError, template.to_netcdf, self.temp_nc_file)
        self.assertIsNone(template.ncobj)
        # self.assertFalse(template.ncobj.isopen())
        # TODO: Use mock to make this fail *after* ncobj is created

    def test_dimensionless_variable(self):
        template = DatasetTemplate(variables={'X': {'_datatype': 'double', '_data': np.array(1)}})
        template.to_netcdf(self.temp_nc_file)

        dataset = Dataset(self.temp_nc_file)
        self.assertEqual((), dataset.variables['X'].dimensions)

    def test_ensure_completeness(self):
        template = DatasetTemplate(dimensions={'X': 1})
        template.variables = {
            'A': {'_dimensions': ['X'], '_datatype': 'float32', '_data': [12.3]},
            'B': {'_dimensions': ['X'], '_data': [12.3]},
            'X': {'_dimensions': ['X'], '_data': self.values1},
            'Y': {'_datatype': 'float32', '_data': None}
        }
        template.ensure_completeness()

        self.assertEqual(['X'], template.variables['A']['_dimensions'])
        self.assertEqual(np.dtype('float32'), template.variables['A']['_datatype'])
        self.assertEqual([12.3], template.variables['A']['_data'])
        self.assertIsInstance(template.variables['A']['_data'], np.ndarray)

        self.assertEqual(np.dtype('float64'), template.variables['B']['_datatype'])

        self.assertIs(self.values1.dtype, template.variables['X']['_datatype'])

        self.assertEqual([], template.variables['Y']['_dimensions'])

        template.variables = {'Z': {'_dimensions': [], '_data': None}}
        self.assertRaisesRegexp(ValidationError, r"No data type information for variable 'Z'",
                                template.ensure_completeness)

        template.variables = {'Z': {'_dimensions': []}}
        self.assertRaisesRegexp(ValidationError, r"No data specified for variable 'Z'",
                                template.ensure_completeness)

    def test_ensure_consistency(self):
        template = DatasetTemplate()
        scalar = {'_dimensions': [], '_data': np.array(1)}
        template.variables = {'SCALAR': scalar}
        template.ensure_consistency()
        self.assertEqual({}, template.dimensions)
        self.assertIs(scalar, template.variables['SCALAR'])

        template = DatasetTemplate(dimensions={'TEN': 10})
        var_10 = {'_dimensions': ['TEN'], '_data': self.values10}
        template.variables = {'TEN': var_10}
        template.ensure_consistency()
        self.assertEqual({'TEN': 10}, template.dimensions)
        self.assertIs(var_10, template.variables['TEN'])

        template = DatasetTemplate(dimensions={'X': None})
        var_12 = {'_dimensions': ['X'], '_data': np.arange(12)}
        template.variables = {'X': var_12}
        template.ensure_consistency()
        self.assertEqual({'X': 12}, template.dimensions)
        self.assertIs(var_12, template.variables['X'])

        empty = {'_dimensions': ['X'], '_data': None}
        template.variables['EMPTY'] = empty
        template.ensure_consistency()
        self.assertEqual({'X': 12}, template.dimensions)
        self.assertIs(empty, template.variables['EMPTY'])

        template.variables['X']['_data'] = self.values1
        self.assertRaisesRegexp(ValueError, 'inconsistent with dimension sizes defined in template',
                                template.ensure_consistency)  # now should fail because dim X is already set

        template.variables = {
            'Z': {'_dimensions': ["NOSUCHTHING"], '_data': self.values10}
        }
        self.assertRaisesRegexp(ValidationError, 'undefined dimensions', template.ensure_consistency)

        template.variables = {
            'W': {'_dimensions': ['X'], '_data': np.arange(20).reshape((10,2))}
        }
        self.assertRaisesRegexp(ValueError,
                                "Variable 'W' has 1 dimensions, but value array has 2 dimensions.",
                                template.ensure_consistency
                                )


class TestDataValues(TemplateTestCase):
    def setUp(self):
        super(TestDataValues, self).setUp()
        self.data_array = np.array([-999., -999., -999., -999., -999., 1., 2., 3., 4., 5])
        self.data_masked = np.ma.masked_array([-4, -3, -2, -1, 0, 1., 2., 3., 4., 5],
                                              mask=[True, True, True, True, True, False, False, False, False, False])
        self.template = DatasetTemplate(
            dimensions={'TIME': 10},
            variables={
                'TIME': {
                    '_dimensions': ['TIME'],
                    '_datatype': 'float64',
                    'valid_min': 0,
                    'valid_max': 10,
                    '_data': np.array([np.nan, np.nan, 1, 2, 3, 4, 5, 6, 7, 8])
                },
                'X': {
                    '_dimensions': ['TIME'],
                    '_datatype': 'float32',
                    'valid_min': 1,
                    'valid_max': 5,
                    '_FillValue': -999,
                    '_data': self.data_array
                },
                'Y': {
                    '_dimensions': ['TIME'],
                    '_datatype': 'float32',
                    'valid_range': [-4, 5],
                    '_fill_value': -999,
                    '_data': self.data_masked
                },
                'N': {
                    '_dimensions': ['TIME'],
                    '_datatype': 'int32',
                    'valid_range': [-4, 5],
                    '_fill_value': -999,
                    '_data': self.data_array
                }
            }
        )

    def test_fill_values(self):
        self.template.to_netcdf(self.temp_nc_file)
        dataset = Dataset(self.temp_nc_file)
        dataset.set_auto_mask(True)
        for varname in ('X', 'Y'):
            dsvar = dataset.variables[varname]
            self.assertEqual(-999., dsvar._FillValue)
            self.assertIsInstance(dsvar[:], np.ma.MaskedArray)
            self.assertTrue(dsvar[:5].mask.all())
            self.assertTrue((dsvar[5:] == self.data_array[5:]).all())

    def test_fill_value_aliases(self):
        self.template.variables['X']['_fill_value'] = -999.  # both aliases, but equal so should still work
        self.template.to_netcdf(self.temp_nc_file)
        dataset = Dataset(self.temp_nc_file)
        self.assertEqual(-999., dataset.variables['X']._FillValue)

        del self._temp_nc_file
        self.template.variables['X']['_fill_value'] = -666.  # now they're different, which is an error
        self.assertRaises(ValueError, self.template.to_netcdf, self.temp_nc_file)

    def test_get_data_range(self):
        self.assertEqual((1, 8), self.template.get_data_range('TIME'))
        self.assertEqual((1, 5), self.template.get_data_range('X'))
        self.assertEqual((1, 5), self.template.get_data_range('Y'))

    def test_var_attr_datatype_conversion(self):
        """
        test to check the conversion of some attributes matches the datatype of the variable as
        defined in the template
        """
        self.template.to_netcdf(self.temp_nc_file)
        dataset = Dataset(self.temp_nc_file)

        TIME = dataset.variables['TIME']
        self.assertEqual(TIME.dtype, TIME.valid_min.dtype)
        self.assertEqual(TIME.dtype, TIME.valid_max.dtype)

        X = dataset.variables['X']
        self.assertEqual(X.dtype, X.valid_min.dtype)
        self.assertEqual(X.dtype, X.valid_max.dtype)
        self.assertEqual(X.dtype, X._FillValue.dtype)

        for v in ['Y', 'N']:
            var = dataset.variables[v]
            self.assertEqual(var.dtype, var.valid_range.dtype)
            self.assertEqual(var.dtype, var._FillValue.dtype)

# TODO: add data from multiple numpy arrays
# e.g. template.add_data(TIME=time_values, TEMP=temp_values, PRES=pres_values)
# TODO: add data from Pandas dataframe (later...)
# e.g. template.add_data(dataframe)

# TODO: create netCDF file with auto-generated file name according to IMOS conventions
# e.g. template.create()


if __name__ == '__main__':
    unittest.main()
