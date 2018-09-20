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

ATTRIBUTES_SCHEMA = {
    "type": "object",
    "patternProperties": {
        NAME_PATTERN: {
            "type": ["string", "number", "array"]
        }
    },
    "additionalProperties": False
}


VARIABLE_DEFINITION_SCHEMA = {
    "type": "object",
    "properties": {
        "dimensions": {
            "type": "array",
            "items": {"type": "string", "pattern": NAME_PATTERN}
        },
        "type": {"type": "string"},
        "attributes": ATTRIBUTES_SCHEMA
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


def validate_dimensions(d):
    validate(d, DIMENSIONS_SCHEMA)


def validate_variables(v):
    validate(v, VARIABLES_SCHEMA)


def validate_attributes(a):
    validate(a, ATTRIBUTES_SCHEMA)
