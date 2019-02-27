from datetime import datetime

from ncwriter import DatasetTemplate


class ImosTemplate(DatasetTemplate):
    """Template class to create IMOS-compliant netCDF files"""

    TIMESTAMP_FORMAT = '%Y-%m-%dT%H:%M:%SZ'

    def __init__(self, *args, **kwargs):
        super(ImosTemplate, self).__init__(*args, **kwargs)
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

