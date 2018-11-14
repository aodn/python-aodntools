"""This module holds schema definitions for validating the various :py:class:`dicts` that make up parts of a
template, and also the helper functions necessary to validate an object against their respective schema.
"""

from jsonschema import validate, ValidationError


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
        "_datatype": {"type": "string"},
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


def validate_dimensions(d):
    validate(d, DIMENSIONS_SCHEMA)


def validate_variables(v):
    validate(v, VARIABLES_SCHEMA)


def validate_global_attributes(a):
    validate(a, ATTRIBUTES_SCHEMA)
