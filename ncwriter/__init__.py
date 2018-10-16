from .template import DatasetTemplate
from .schema import ValidationError

__all__ = [
    'DatasetTemplate',
    'ValidationError'
]

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
