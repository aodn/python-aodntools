#!/usr/bin/env python3

import os
import unittest

import xarray as xr

from aodntools.timeseries_products.common import (check_file, check_velocity_file, get_qc_variable_names,
                                                  check_imos_flag_conventions, in_water_index, in_water)

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
AM_FILE = os.path.join(
    TEST_ROOT, 'IMOS_ANMN-AM_GST_20190419T100000Z_NRSMAI_FV01_NRSMAI-CO2-1904-delayed_END-20190531T020000Z_C-20200625T000000Z.nc'
)


class TestQCVariableFunctions(unittest.TestCase):
    def test_get_qc_variable_names(self):
        with xr.open_dataset(GOOD_TZ_FILE) as nc:
            qc_var = get_qc_variable_names(nc)
        self.assertEqual(qc_var, ['DEPTH_quality_control', 'TEMP_quality_control'])

    def test_check_flag_convetions(self):
        with xr.open_dataset(GOOD_TZ_FILE) as nc:
            self.assertEqual(check_imos_flag_conventions(nc), [])
            self.assertEqual(check_imos_flag_conventions(nc, ['NONE']),
                             ['variable NONE not in file'])
            self.assertEqual(check_imos_flag_conventions(nc, ['TEMP']),
                             ['variable TEMP missing quality_control_conventions'])

    def test_check_flag_conventions_bad(self):
        with xr.open_dataset(AM_FILE) as nc:
            errors = check_imos_flag_conventions(nc)
        self.assertEqual(errors, ['unexpected quality_control_conventions: "WOCE quality control procedure"'])


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
            error_list = check_velocity_file(nc, 'NRSROT')
        self.assertEqual(error_list, [])

    def test_bad_velocity_file(self):
        with xr.open_dataset(BAD_V_FILE) as nc:
            error_list = check_velocity_file(nc, 'NWSROW')
        self.assertEqual(set(error_list), {'VCUR variable missing',
                                           'DEPTH variable missing',
                                           "dimension(s) {'DIST_ALONG_BEAMS'} not allowed for UCUR",
                                           'no in-water data'
                                           }
                         )

    def test_am_file(self):
        with xr.open_dataset(AM_FILE) as nc:
            error_list = check_file(nc, 'NRSMAI', 'TEMP')
        self.assertEqual(set(error_list), {'no NOMINAL_DEPTH',
                                           'no time_deployment_start attribute',
                                           'no time_deployment_end attribute',
                                           'unexpected quality_control_conventions: "WOCE quality control procedure"'
                                           }
                         )


class TestInWater(unittest.TestCase):
    def test_in_water_index_ok(self):
        with xr.open_dataset(BAD_TZ_FILE) as nc:
            nc.attrs['time_deployment_start'] = '2018-12-13T08:00:00Z'
            nc.attrs['time_deployment_end'] = '2018-12-14T00:30:00Z'
            index = in_water_index(nc)
        self.assertTrue(all(index[:-2]))
        self.assertFalse(any(index[-2:]))

    def test_in_water_index_bad(self):
        with xr.open_dataset(BAD_V_FILE) as nc:
            index = in_water_index(nc)
        self.assertFalse(all(index))

    def test_in_water_ok(self):
        with xr.open_dataset(BAD_TZ_FILE) as nc:
            nc.attrs['time_deployment_start'] = '2018-12-13T08:00:00Z'
            nc.attrs['time_deployment_end'] = '2018-12-14T00:30:00Z'
            nc_in = in_water(nc)

            self.assertEqual(len(nc_in.TIME), len(nc.TIME) - 2)
            self.assertTrue(all(nc_in.TIME.values == nc.TIME[:-2].values))

if __name__ == '__main__':
    unittest.main()
