#!/usr/bin/env python3

import unittest
from datetime import datetime, timedelta

import numpy as np

from aodntools.ncwriter import ImosTemplate
from aodntools.ncwriter import DatasetTemplate
from test_aodntools.ncwriter.test_template import TemplateTestCase

TEST_FIXED_GLOBALS = {
    "Conventions": "CF-1.6,IMOS-1.4",
    "project": "Integrated Marine Observing System (IMOS)",
    "naming_authority": "IMOS",
    "data_centre": "Australian Ocean Data Network (AODN)",
    "data_centre_email": "info@aodn.org.au"
}


class TestImosTemplate(TemplateTestCase):

    def setUp(self):
        super(TestImosTemplate, self).setUp()
        self.template = ImosTemplate()

    def test_init(self):
        self.assertIsInstance(self.template, DatasetTemplate)

    def test_fixed_global_attributes(self):
        for name, value in TEST_FIXED_GLOBALS.items():
            self.assertEqual(value, self.template.global_attributes[name])

    def test_combine_global_attributes(self):
        my_globals = {"title": "This is a test",
                      "project": "Test project"
                      }
        expected = TEST_FIXED_GLOBALS.copy()
        expected.update(my_globals)
        template = ImosTemplate(global_attributes=my_globals)

        for name, value in expected.items():
            self.assertEqual(value, template.global_attributes[name])

    def test_date_created(self):
        now = datetime.now()
        self.assertIsInstance(self.template.date_created, datetime)
        self.assertTrue(self.template.date_created - now < timedelta(seconds=1))

    def test_add_date_created_attribute(self):
        self.template.add_date_created_attribute()
        parsed_date_created = datetime.strptime(self.template.global_attributes['date_created'], '%Y-%m-%dT%H:%M:%SZ')
        self.assertTrue(self.template.date_created - parsed_date_created < timedelta(seconds=1))

    def test_add_extent_attributes(self):
        self.template.add_extent_attributes(time_var=None, vert_var=None, lat_var=None, lon_var=None)
        for att in ('time_coverage_start', 'time_coverage_end', 'geospatial_vertical_min', 'geospatial_vertical_max',
                    'geospatial_lat_min', 'geospatial_lat_max', 'geospatial_lon_min', 'geospatial_lon_max'):
            self.assertEqual('', self.template.global_attributes[att])

        self.template.variables = {
            'TIME': {
                '_dimensions': ['TIME'],
                '_datatype': 'float64',
                '_data': np.array([0, 1, 2, np.nan, 25261.375]),
                'units': 'days since 1950-01-01 00:00:00 UTC'
            },
            'LATITUDE': {
                '_dimensions': ['TIME'],
                '_datatype': 'float32',
                '_FillValue': -999.,
                '_data': np.array([-999., -999., -42, -43, 12])
            },
            'LONGITUDE': {
                '_dimensions': ['TIME'],
                '_datatype': 'float32',
                '_data': np.arange(10)
            },
            'NOMINAL_DEPTH': {
                '_datatype': 'float32',
                '_data': 20
            },
            'DEPTH': {
                '_dimensions': ['TIME'],
                '_datatype': 'float32',
                '_data': np.repeat(np.nan, 5)
            }
        }
        self.template.add_extent_attributes(vert_var='NOMINAL_DEPTH')
        self.assertEqual('1950-01-01T00:00:00Z', self.template.global_attributes['time_coverage_start'])
        self.assertEqual('2019-03-01T09:00:00Z', self.template.global_attributes['time_coverage_end'])
        self.assertEqual(-43, self.template.global_attributes['geospatial_lat_min'])
        self.assertEqual(12, self.template.global_attributes['geospatial_lat_max'])
        self.assertEqual(0, self.template.global_attributes['geospatial_lon_min'])
        self.assertEqual(9, self.template.global_attributes['geospatial_lon_max'])
        self.assertEqual(20, self.template.global_attributes['geospatial_vertical_min'])
        self.assertEqual(20, self.template.global_attributes['geospatial_vertical_max'])

        self.assertRaises(ValueError, self.template.add_extent_attributes)


if __name__ == '__main__':
    unittest.main()
