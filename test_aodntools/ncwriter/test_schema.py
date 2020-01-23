#!/usr/bin/env python3

import unittest

import numpy as np

from aodntools.ncwriter.schema import (validate_template, validate_dimensions, validate_variables,
                                       validate_global_attributes, ValidationError)


class TestSchema(unittest.TestCase):
    def test_validate_template(self):
        validate_template({})
        validate_template({'_dimensions': {'X': 1}})
        validate_template({'_variables': {'X': {'name': 'X'}}})
        validate_template({'title': 'test'})

        self.assertRaises(ValidationError, validate_template, {'_bad': 1})
        self.assertRaises(ValidationError, validate_template, {'_dimensions': None})
        self.assertRaises(ValidationError, validate_template, {'_variables': 1})
        self.assertRaises(ValidationError, validate_template, {'_variables': {'X': {'_dimensions': 1}}})

    def test_validate_dimensions(self):
        validate_dimensions({})
        validate_dimensions({'X': 1, 'Y': None})

        with self.assertRaises(ValidationError):
            validate_dimensions(None)
        with self.assertRaises(ValidationError):
            validate_dimensions('X')
        with self.assertRaises(ValidationError):
            validate_dimensions(['X'])
        with self.assertRaises(ValidationError):
            validate_dimensions(10)

        with self.assertRaises(ValidationError):
            validate_dimensions({'123': 123})
        with self.assertRaises(ValidationError):
            validate_dimensions({'X': 'one'})
        with self.assertRaises(ValidationError):
            validate_dimensions({'X': -1})
        with self.assertRaises(ValidationError):
            validate_dimensions({'X': 1.5})

    def test_validate_variables(self):
        validate_variables({})
        validate_variables({'X': {}})
        validate_variables({'X': {'_datatype': 'float32'}})
        validate_variables({'X': {'_datatype': np.float32}})
        validate_variables({'X': {'_datatype': np.dtype('float32')}})
        validate_variables({'X': {'_dimensions': []}})
        validate_variables({'X': {'name': 'X'}})
        validate_variables({'X': {'_data': None}})
        validate_variables({'X': {'_data': 100}})
        validate_variables({'X': {'_data': np.array([1, 2])}})
        validate_variables({'X': {'_dimensions': ['X'], '_datatype': 'float32'}})
        validate_variables({
            'X': {
                '_dimensions': ['X'],
                '_datatype': 'float32',
                '_data': [42],
                'name': 'X',
                'count': 1
            },
            'Y': {'_datatype': 'float64'}
        })

        with self.assertRaises(ValidationError):
            validate_variables(None)
        with self.assertRaises(ValidationError):
            validate_variables('VAR')
        with self.assertRaises(ValidationError):
            validate_variables({'__X': {}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'_unknown': 'else'}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'_datatype': 'no_such_type'}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'_datatype': 42}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'_datatype': 'float32', '0': 'none'}})

    def test_validate_attributes(self):
        validate_global_attributes({})
        validate_global_attributes({'name': 'test'})
        validate_global_attributes({'name': 1.5})
        validate_global_attributes({'name': [1, 2, 3]})

        with self.assertRaises(ValidationError):
            validate_global_attributes(None)
        with self.assertRaises(ValidationError):
            validate_global_attributes('X')
        with self.assertRaises(ValidationError):
            validate_global_attributes([])
        with self.assertRaises(ValidationError):
            validate_global_attributes({'_dimensions': {}})
        with self.assertRaises(ValidationError):
            validate_global_attributes({'_badname': 1})
        with self.assertRaises(ValidationError):
            validate_global_attributes({'null': None})
        with self.assertRaises(ValidationError):
            validate_global_attributes({'bool': True})
        with self.assertRaises(ValidationError):
            validate_global_attributes({'object': {}})


if __name__ == '__main__':
    unittest.main()
