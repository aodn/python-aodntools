"""This module holds schema definitions for validating the various :py:class:`dicts` that make up parts of a
template, and also the helper functions necessary to validate an object against their respective schema.
"""

import json

import numpy as np

from jsonschema import validators, Draft4Validator, FormatChecker, ValidationError
from pkg_resources import resource_filename

try:
    basestring
except NameError:
    basestring = str

# Create a new validator class (based on Draft4Validator) to allow templates to use
# * Python types or numpy dtypes to specify variable data types; and
# * numpy arrays to specify variable data.

TemplateValidator = validators.create(meta_schema=Draft4Validator.META_SCHEMA,
                                      validators=Draft4Validator.VALIDATORS)
format_checker = FormatChecker()

TYPE_STRINGS = ('int', 'int8', 'int16', 'int32', 'int64',
                'float', 'float8', 'float16', 'float32', 'float64', 'double',
                'S1', 'c', 'i1', 'b', 'B', 'u1', 'i2', 'h', 's', 'u2', 'i4', 'i', 'l', 'u4',
                'i8', 'u8', 'f4', 'f', 'f8', 'd'
                )


@format_checker.checks('datatype')
def is_datatype(value):
    """Return whether the given value is a valid data type specification for a NetCDF variable"""
    if isinstance(value, basestring):
        return value in TYPE_STRINGS
    if isinstance(value, np.dtype):
        return True
    if isinstance(value, type):
        return issubclass(value, np.number)

    return False


TYPES = {'array': (list, np.ndarray)}

TEMPLATE_SCHEMA_JSON = resource_filename(__name__, 'template_schema.json')
with open(TEMPLATE_SCHEMA_JSON) as f:
    TEMPLATE_SCHEMA = json.load(f)
TemplateValidator.check_schema(TEMPLATE_SCHEMA)

template_validator = TemplateValidator(TEMPLATE_SCHEMA, types=TYPES, format_checker=format_checker)

GLOBAL_ATTRIBUTES_SCHEMA = TEMPLATE_SCHEMA.copy()
GLOBAL_ATTRIBUTES_SCHEMA.pop('properties')  # remove special properties, leaving only global attributes
global_attribute_validator = TemplateValidator(GLOBAL_ATTRIBUTES_SCHEMA, types=TYPES, format_checker=format_checker)


def validate_template(t):
    template_validator.validate(t)


def validate_dimensions(d):
    validate_template({'_dimensions': d})


def validate_variables(v):
    validate_template({'_variables': v})


def validate_global_attributes(a):
    global_attribute_validator.validate(a)
