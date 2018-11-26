"""This module holds schema definitions for validating the various :py:class:`dicts` that make up parts of a
template, and also the helper functions necessary to validate an object against their respective schema.
"""

from copy import deepcopy

import numpy as np

from jsonschema import validators, Draft4Validator, ValidationError

try:
    basestring
except NameError:
    basestring = str

# Create a new validator class (based on Draft4Validator) to allow templates to use
# * Python types or numpy dtypes to specify variable data types; and
# * numpy arrays to specify variable data.
meta = deepcopy(Draft4Validator.META_SCHEMA)
meta['definitions']['simpleTypes']['enum'].append('datatype')
types = deepcopy(Draft4Validator.DEFAULT_TYPES)
types['datatype'] = (basestring, type, np.dtype)
types['array'] = (list, np.ndarray)
TemplateValidator = validators.create(meta_schema=meta,
                                      validators=Draft4Validator.VALIDATORS,
                                      default_types=types
                                      )


NAME_PATTERN = r"^[A-Za-z][A-Za-z0-9_]*$"


DIMENSIONS_SCHEMA = {
    "type": "object",
    "patternProperties": {
        NAME_PATTERN: {
            "type": ["integer", "null"],
            "minimum": 0
        }
    },
    "additionalProperties": False
}

VARIABLE_DEFINITION_SCHEMA = {
    "type": "object",
    "properties": {
        "_dimensions": {
            "type": "array",
            "items": {"type": "string", "pattern": NAME_PATTERN}
        },
        "_datatype": {"type": "datatype"},
        "_FillValue": {"type": ["number", "string"]},
        "_data": {"type": ["null", "array"]}
    },
    "patternProperties": {
        NAME_PATTERN: {"type": ["string", "number", "array"]}
    },
    "additionalProperties": False
}


VARIABLES_SCHEMA = {
    "type": "object",
    "patternProperties": {
        NAME_PATTERN: VARIABLE_DEFINITION_SCHEMA
    },
    "additionalProperties": False
}

TEMPLATE_SCHEMA = {
    "type": "object",
    "properties": {
        "_dimensions": DIMENSIONS_SCHEMA,
        "_variables": VARIABLES_SCHEMA
    },
    "patternProperties": {
        NAME_PATTERN: {
            "type": ["string", "number", "array"]
        }
    },
    "additionalProperties": False
}


ATTRIBUTES_SCHEMA = TEMPLATE_SCHEMA.copy()
ATTRIBUTES_SCHEMA.pop("properties")  # remove special properties, leaving only global attributes


TemplateValidator.check_schema(TEMPLATE_SCHEMA)


def validate_dimensions(d):
    TemplateValidator(DIMENSIONS_SCHEMA).validate(d)


# validate_variables = LocalValidator(VARIABLES_SCHEMA).validate
def validate_variables(v):
    TemplateValidator(VARIABLES_SCHEMA).validate(v)


def validate_global_attributes(a):
    TemplateValidator(ATTRIBUTES_SCHEMA).validate(a)
