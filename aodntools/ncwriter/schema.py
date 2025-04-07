"""This module holds schema definitions for validating the various :py:class:`dicts` that make up parts of a
template, and also the helper functions necessary to validate an object against their respective schema.
"""
import json
import numpy as np
from jsonschema import validators, Draft4Validator, FormatChecker, ValidationError
from pkg_resources import resource_filename

def is_array(checker, instance):
    return isinstance(instance, (list, np.ndarray))

# Extend the default type checker by redefining "array"
custom_type_checker = Draft4Validator.TYPE_CHECKER.redefine("array", is_array)

# Create a custom validator that uses the new type checker.
CustomValidator = validators.extend(Draft4Validator, type_checker=custom_type_checker)
format_checker = FormatChecker()

@format_checker.checks('datatype')
def is_python_datatype(value):
    """Return whether the given value is a valid data type specification for a NetCDF variable"""
    if isinstance(value, np.dtype):
        return True
    if isinstance(value, type):
        return issubclass(value, np.number)
    return False

TEMPLATE_SCHEMA_JSON = resource_filename(__name__, 'template_schema.json')
with open(TEMPLATE_SCHEMA_JSON) as f:
    TEMPLATE_SCHEMA = json.load(f)
CustomValidator.check_schema(TEMPLATE_SCHEMA)

# Use the custom validator
template_validator = CustomValidator(TEMPLATE_SCHEMA, format_checker=format_checker)

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
