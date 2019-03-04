from datetime import datetime
from pkg_resources import resource_filename

from ncwriter import DatasetTemplate


IMOS_GLOBAL_JSON = resource_filename(__name__, 'imos_global.json')
IMOS_GLOBAL_ATTRIBUTES = DatasetTemplate.from_json(IMOS_GLOBAL_JSON).global_attributes


class ImosTemplate(DatasetTemplate):
    """Template class to create IMOS-compliant netCDF files"""

    TIMESTAMP_FORMAT = '%Y-%m-%dT%H:%M:%SZ'

    def __init__(self, global_attributes=None, *args, **kwargs):
        # add standard IMOS global attributes first
        combined_attributes = IMOS_GLOBAL_ATTRIBUTES.copy()
        if global_attributes is not None:
            combined_attributes.update(global_attributes)
        super(ImosTemplate, self).__init__(global_attributes=combined_attributes, *args, **kwargs)
        self._date_created = None

    @property
    def date_created(self):
        """Read-only property to return the creation date of the netCDF file as a datetime object"""
        if self._date_created is None:
            self._date_created = datetime.now()
        return self._date_created

    def add_date_created_attribute(self):
        self.global_attributes['date_created'] = self.date_created.strftime(self.TIMESTAMP_FORMAT)

    # TODO: def add_extent_attributes(self):

    # TODO: def set_imos_filename(self):

