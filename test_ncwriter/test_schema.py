import unittest

from ncwriter.schema import validate_dimensions, validate_variables, validate_global_attributes, ValidationError


class TestSchema(unittest.TestCase):
    def test_validate_dimensions(self):
        validate_dimensions({})
        validate_dimensions({'X': 1, 'Y': 2})

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
        validate_variables({'X': {'_dimensions': []}})
        validate_variables({'X': {'name': 'X'}})
        validate_variables({'X': {'_data': None}})
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
            validate_variables({'X': {'_datatype': 'float32', '_unknown': 'else'}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'_datatype': 'float32', '0': 'none'}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'_data': 100}})

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
            validate_global_attributes({'_badname': 1})
        with self.assertRaises(ValidationError):
            validate_global_attributes({'null': None})
        with self.assertRaises(ValidationError):
            validate_global_attributes({'bool': True})
        with self.assertRaises(ValidationError):
            validate_global_attributes({'object': {}})


if __name__ == '__main__':
    unittest.main()
