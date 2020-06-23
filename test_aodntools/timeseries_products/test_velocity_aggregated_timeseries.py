#!/usr/bin/env python3

import os
import unittest

from netCDF4 import Dataset, chartostring

from aodntools import __version__
from aodntools.timeseries_products.velocity_aggregated_timeseries import velocity_aggregated
from test_aodntools.base_test import BaseTestCase

TEST_ROOT = os.path.dirname(__file__)
INPUT_FILES = [
    'IMOS_ANMN-NRS_AETVZ_20181213T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1812-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20181215T100000Z_C-20200430T000000Z.nc',
    'IMOS_ANMN-NRS_AETVZ_20180816T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1808-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20180822T053000Z_C-20200623T000000Z.nc',
    'IMOS_ANMN-NRS_AETVZ_20191016T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1910-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20191018T100000Z_C-20200430T000000Z.nc',
]
EXPECTED_OUTPUT_FILE = os.path.join(
    TEST_ROOT, 'IMOS_ANMN-NRS_VZ_20180816_NRSROT_FV01_velocity-aggregated-timeseries_END-20191018_C-20200623.nc'
)

class TestVelocityAggregatedTimeseries(BaseTestCase):
    def test_main_aggregator(self):
        output_file, bad_files = velocity_aggregated(INPUT_FILES, 'NRSROT', input_dir=TEST_ROOT, output_dir='/tmp')

        self.assertEqual(0, len(bad_files))

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


if __name__ == '__main__':
    unittest.main()
