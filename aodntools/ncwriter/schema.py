"""This module holds schema definitions for validating the various :py:class:`dicts` that make up parts of a
template, and also the helper functions necessary to validate an object against their respective schema.
"""
import json
import numpy as np
from jsonschema import validators, Draft4Validator, FormatChecker, ValidationError
from pathlib import Path

# helper function that will later be used to tell the schema validator how to validate objects of type "array"
def is_array(checker, instance):
    return isinstance(instance, (list, np.ndarray))

# Extend the default type checker by redefining "array"
# whenever a schema expects a value of type "array", it will now use the is_array function to check if the value is acceptable.
custom_type_checker = Draft4Validator.TYPE_CHECKER.redefine("array", is_array)

# Create a custom validator that uses the new type checker.
# any validation performed with CustomValidator will use the custom array checker 
CustomValidator = validators.extend(Draft4Validator, type_checker=custom_type_checker)
format_checker = FormatChecker()

# Define a custom format checker
# called when a JSON schema specifies that a value should have the format "datatype"
@format_checker.checks('datatype')
def is_python_datatype(value):
    """Return whether the given value is a valid data type specification for a NetCDF variable"""
    if isinstance(value, np.dtype):
        return True
    if isinstance(value, type):
        return issubclass(value, np.number)
    return False

# Load JSON schema file
TEMPLATE_SCHEMA_JSON = Path(__file__).parent  / 'template_schema.json'
with open(TEMPLATE_SCHEMA_JSON) as f:
    TEMPLATE_SCHEMA = json.load(f)

# Use the custom validator to check it is valid according to Draft 4 rules
CustomValidator.check_schema(TEMPLATE_SCHEMA)

# ready-to-use validator that applies both custom type and format checks
template_validator = CustomValidator(TEMPLATE_SCHEMA, format_checker=format_checker)


# Validation checks
def validate_template(t):
    template_validator.validate(t)

def validate_dimensions(d):
    validate_template({'_dimensions': d})

def validate_variables(v):
    validate_template({'_variables': v})
    
def validate_global_attributes(a):
    if hasattr(a, 'keys'):
        special = [k for k in a.keys() if k.startswith('_')]
        if special:
            raise ValidationError('Special attributes {} not allowed in global attributes dict'.format(special))
    template_validator.validate(a)
