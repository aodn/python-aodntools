#!/usr/bin/env python3

import os
import unittest

from netCDF4 import Dataset, chartostring

from aodntools import __version__
from aodntools.timeseries_products.common import NoInputFilesError
from aodntools.timeseries_products.velocity_aggregated_timeseries import velocity_aggregated
from test_aodntools.base_test import BaseTestCase

TEST_ROOT = os.path.dirname(__file__)
BAD_FILE = 'IMOS_ANMN-NRS_TZ_20181213T080000Z_NRSROT_FV00_NRSROT-1812-SBE39-43_END-20181214T004000Z_C-20190827T000000Z.nc'
INPUT_FILES = [
    'IMOS_ANMN-NRS_AETVZ_20181213T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1812-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20181215T100000Z_C-20200430T000000Z.nc',
    'IMOS_ANMN-NRS_AETVZ_20180816T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1808-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20180822T053000Z_C-20200623T000000Z.nc',
    'IMOS_ANMN-NRS_AETVZ_20191016T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1910-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20191018T100000Z_C-20200430T000000Z.nc',
    BAD_FILE
]
EXPECTED_OUTPUT_FILE = os.path.join(
    TEST_ROOT, 'IMOS_ANMN-NRS_VZ_20180816_NRSROT_FV01_velocity-aggregated-timeseries_END-20191018_C-20200623.nc'
)

class TestVelocityAggregatedTimeseries(BaseTestCase):
    def test_main_aggregator(self):
        output_file, bad_files = velocity_aggregated(INPUT_FILES, 'NRSROT', input_dir=TEST_ROOT, output_dir='/tmp')

        self.assertEqual(1, len(bad_files))
        for file, errors in bad_files.items():
            self.assertEqual(BAD_FILE, file)
            self.assertSetEqual(set(errors), {'no NOMINAL_DEPTH',
                                              'Wrong file version: Level 0 - Raw Data',
                                              'DEPTH variable missing',
                                              'UCUR variable missing',
                                              'VCUR variable missing',
                                              'no time_deployment_start attribute',
                                              'no time_deployment_end attribute'
                                              }
                                )

        dataset = Dataset(output_file)

        # check dimensions and variables
        self.assertSetEqual(set(dataset.dimensions), {'OBSERVATION', 'INSTRUMENT', 'strlen'})

        obs_vars = {'TIME', 'DEPTH', 'DEPTH_quality_control', 'UCUR', 'UCUR_quality_control',
                    'VCUR', 'VCUR_quality_control', 'WCUR', 'WCUR_quality_control', 'instrument_index', 'CELL_INDEX'}
        inst_vars = {'LATITUDE', 'LONGITUDE', 'NOMINAL_DEPTH', 'SECONDS_TO_MIDDLE'}
        string_vars = {'source_file', 'instrument_id'}

        self.assertSetEqual(set(dataset.variables.keys()), obs_vars | inst_vars | string_vars)

        for var in obs_vars:
            self.assertEqual(dataset.variables[var].dimensions, ('OBSERVATION',))
        for var in inst_vars:
            self.assertEqual(dataset.variables[var].dimensions, ('INSTRUMENT',))
        for var in string_vars:
            self.assertEqual(dataset.variables[var].dimensions, ('INSTRUMENT', 'strlen'))

        for f in chartostring(dataset['source_file'][:]):
            self.assertIn(f, INPUT_FILES)

        # check attributes
        self.assertEqual(__version__, dataset.generating_code_version)
        self.assertIn(__version__, dataset.lineage)

        # check aggregated variable values
        expected = Dataset(EXPECTED_OUTPUT_FILE)
        non_match_vars = []
        for var in ('TIME', 'UCUR', 'UCUR_quality_control', 'VCUR', 'VCUR_quality_control', 'NOMINAL_DEPTH',
                    'CELL_INDEX', 'instrument_index', 'DEPTH', 'DEPTH_quality_control'):
            self.assertTrue(all(dataset[var][:] == expected[var][:]), "{} values don't match up!".format(var))

    def test_all_rejected(self):
        self.assertRaises(NoInputFilesError, velocity_aggregated, [BAD_FILE], 'NRSROT', input_dir=TEST_ROOT)


if __name__ == '__main__':
    unittest.main()
