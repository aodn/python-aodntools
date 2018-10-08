import unittest

from ncwriter.schema import validate_dimensions, validate_variables, validate_attributes, ValidationError


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
        validate_variables({'X': {'type': 'float32'}})
        validate_variables({'X': {'dimensions': []}})
        validate_variables({'X': {'attributes': {'name': 'X'}}})
        validate_variables({'X': {'data': None}})
        validate_variables({'X': {'dimensions': ['X'], 'type': 'float32'}})
        validate_variables({
            'X': {
                'dimensions': ['X'],
                'type': 'float32',
                'attributes': {'name': 'X', 'count': 1},
                'data': [42]
            }
        })

        with self.assertRaises(ValidationError):
            validate_variables(None)
        with self.assertRaises(ValidationError):
            validate_variables('VAR')
        with self.assertRaises(ValidationError):
            validate_variables({'__X': {}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'type': 'float32', 'something': 'else'}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'type': 'float32', 'attributes': 'none'}})
        with self.assertRaises(ValidationError):
            validate_variables({'X': {'data': 100}})

    def test_validate_attributes(self):
        validate_attributes({})
        validate_attributes({'name': 'test'})
        validate_attributes({'name': 1.5})
        validate_attributes({'name': [1, 2, 3]})

        with self.assertRaises(ValidationError):
            validate_attributes(None)
        with self.assertRaises(ValidationError):
            validate_attributes('X')
        with self.assertRaises(ValidationError):
            validate_attributes([])
        with self.assertRaises(ValidationError):
            validate_attributes({'_badname': 1})
        with self.assertRaises(ValidationError):
            validate_attributes({'null': None})
        with self.assertRaises(ValidationError):
            validate_attributes({'bool': True})
        with self.assertRaises(ValidationError):
            validate_attributes({'object': {}})


if __name__ == '__main__':
    unittest.main()
