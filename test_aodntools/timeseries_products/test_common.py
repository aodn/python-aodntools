#!/usr/bin/env python3

import os
import unittest

import xarray as xr

from aodntools.timeseries_products.common import check_file

TEST_ROOT = os.path.dirname(__file__)
GOOD_TZ_FILE = os.path.join(
    TEST_ROOT, 'IMOS_ANMN-NRS_TZ_20181213T080000Z_NRSROT_FV01_NRSROT-1812-SBE39-23_END-20190306T160000Z_C-20190827T000000Z.nc'
)
GOOD_V_FILE = os.path.join(
    TEST_ROOT, 'IMOS_ANMN-NRS_AETVZ_20181213T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1812-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20181215T100000Z_C-20200430T000000Z.nc'
)
BAD_TZ_FILE = os.path.join(
    TEST_ROOT, 'IMOS_ANMN-NRS_TZ_20181213T080000Z_NRSROT_FV00_NRSROT-1812-SBE39-43_END-20181214T004000Z_C-20190827T000000Z.nc'
)
BAD_V_FILE = os.path.join(
    TEST_ROOT, 'IMOS_ANMN-NRS_BAD_VELOCITY_FILE.nc'
)

class TestCheckFile(unittest.TestCase):
    def test_good_temp_file(self):
        with xr.open_dataset(GOOD_TZ_FILE) as nc:
            error_list = check_file(nc, 'NRSROT', 'TEMP')
        self.assertEqual(error_list, [])

    def test_variable_list(self):
        with xr.open_dataset(GOOD_TZ_FILE) as nc:
            error_list = check_file(nc, 'NRSROT', ['TEMP', 'PSAL', 'DEPTH'])
        self.assertEqual(error_list, [])

    def test_wrong_site_and_var(self):
        with xr.open_dataset(GOOD_TZ_FILE) as nc:
            error_list = check_file(nc, 'NO_SITE', 'OTHER')
        self.assertEqual(set(error_list), {'Wrong site_code: NRSROT', 'no variables to aggregate'})

    def test_bad_temp_file(self):
        with xr.open_dataset(BAD_TZ_FILE) as nc:
            error_list = check_file(nc, 'NRSROT', 'TEMP')
        self.assertEqual(set(error_list),
                         {'no NOMINAL_DEPTH', 'Wrong file version: Level 0 - Raw Data',
                         'no time_deployment_start attribute', 'no time_deployment_end attribute'}
                         )

    def test_good_velocity_file(self):
        with xr.open_dataset(GOOD_V_FILE) as nc:
            error_list = check_file(nc, 'NRSROT', 'velocity')
        self.assertEqual(error_list, [])

    def test_bad_velocity_file(self):
        with xr.open_dataset(BAD_V_FILE) as nc:
            error_list = check_file(nc, 'NWSROW', 'velocity')
        self.assertEqual(set(error_list), {'VCUR variable missing',
                                           'DEPTH variable missing',
                                           "dimension(s) {'DIST_ALONG_BEAMS'} not allowed for UCUR"
                                           }
                         )


if __name__ == '__main__':
    unittest.main()
