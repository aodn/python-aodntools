"""This module holds schema definitions for validating the various :py:class:`dicts` that make up parts of a
template, and also the helper functions necessary to validate an object against their respective schema.
"""

import json
import os
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

TEMPLATE_SCHEMA_JSON = os.path.join(os.path.dirname(__file__), 'template_schema.json')
with open(TEMPLATE_SCHEMA_JSON) as f:
    TEMPLATE_SCHEMA = json.load(f)
TemplateValidator.check_schema(TEMPLATE_SCHEMA)


DIMENSIONS_SCHEMA = TEMPLATE_SCHEMA['definitions']['fileDimensions']

_TEMP_DEFS_COPY = deepcopy(TEMPLATE_SCHEMA['definitions'])
VARIABLES_SCHEMA = _TEMP_DEFS_COPY.pop('fileVariables')
VARIABLES_SCHEMA['definitions'] = _TEMP_DEFS_COPY

GLOBAL_ATTRIBUTES_SCHEMA = TEMPLATE_SCHEMA.copy()
GLOBAL_ATTRIBUTES_SCHEMA.pop("properties")  # remove special properties, leaving only global attributes


def validate_template(t):
    TemplateValidator(TEMPLATE_SCHEMA).validate(t)


def validate_dimensions(d):
    TemplateValidator(DIMENSIONS_SCHEMA).validate(d)


def validate_variables(v):
    TemplateValidator(VARIABLES_SCHEMA).validate(v)


def validate_global_attributes(a):
    TemplateValidator(GLOBAL_ATTRIBUTES_SCHEMA).validate(a)
