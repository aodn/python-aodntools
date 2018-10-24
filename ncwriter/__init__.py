from .template import DatasetTemplate, metadata_attributes
from .schema import ValidationError

__all__ = [
    'DatasetTemplate',
    'ValidationError',
    'metadata_attributes'
]

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
