from datetime import datetime
from pkg_resources import resource_filename

from netCDF4 import num2date

from .template import DatasetTemplate

IMOS_GLOBAL_JSON = resource_filename(__name__, 'imos_global.json')
IMOS_GLOBAL_ATTRIBUTES = DatasetTemplate.from_json(IMOS_GLOBAL_JSON).global_attributes

TIMESTAMP_FORMAT = '%Y-%m-%dT%H:%M:%SZ'


class ImosTemplate(DatasetTemplate):
    """Template class to create IMOS-compliant netCDF files"""

    def __init__(self, global_attributes=None, *args, **kwargs):
        # add standard IMOS global attributes first
        combined_attributes = IMOS_GLOBAL_ATTRIBUTES.copy()
        if global_attributes is not None:
            combined_attributes.update(global_attributes)
        super(ImosTemplate, self).__init__(global_attributes=combined_attributes, *args, **kwargs)
        self._date_created = None

    @property
    def date_created(self):
        """Read-only property to return the UTC creation time of the netCDF file as a datetime object"""
        if self._date_created is None:
            self._date_created = datetime.utcnow()
        return self._date_created

    def add_date_created_attribute(self):
        self.global_attributes['date_created'] = self.date_created.strftime(TIMESTAMP_FORMAT)

    def add_extent_attributes(self, time_var='TIME', vert_var='DEPTH', lat_var='LATITUDE', lon_var='LONGITUDE'):
        """
        Calculate spatial and temporal extents from coordinate variables in the template and add/update
        the relevant global attributes. Set an input variable name to None to skip that coordinate.

        :param time_var: Name of the time variable to use (default 'TIME')
        :param vert_var: Name of the vertical coordinate variable (default 'DEPTH')
        :param lat_var: Name of the latitude variable (default 'LATITUDE')
        :param lon_var: Name of the longitude variable (default 'LONGITUDE')
        :return: None
        """
        if time_var:
            data_range = self.get_data_range(time_var)
            time = self.variables[time_var]
            units = time.get('units')
            if not units:
                raise ValueError("Time variable '{time_var}' has no units".format(time_var=time_var))
            calendar = time.get('calendar', 'gregorian')
            time_range = num2date(data_range, units, calendar)
            self.global_attributes['time_coverage_start'] = time_range[0].strftime(TIMESTAMP_FORMAT)
            self.global_attributes['time_coverage_end'] = time_range[1].strftime(TIMESTAMP_FORMAT)

        if vert_var:
            vmin, vmax = self.get_data_range(vert_var)
            self.global_attributes['geospatial_vertical_min'] = vmin
            self.global_attributes['geospatial_vertical_max'] = vmax

        if lat_var:
            vmin, vmax = self.get_data_range(lat_var)
            self.global_attributes['geospatial_lat_min'] = vmin
            self.global_attributes['geospatial_lat_max'] = vmax

        if lon_var:
            vmin, vmax = self.get_data_range(lon_var)
            self.global_attributes['geospatial_lon_min'] = vmin
            self.global_attributes['geospatial_lon_max'] = vmax

    # TODO: def set_imos_filename(self):

