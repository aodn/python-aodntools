#!/usr/bin/env python3

import os
import unittest

from netCDF4 import Dataset, chartostring

from aodntools import __version__
from aodntools.timeseries_products.common import NoInputFilesError
from aodntools.timeseries_products.velocity_aggregated_timeseries import velocity_aggregated
from test_aodntools.base_test import BaseTestCase

TEST_ROOT = os.path.dirname(__file__)
BAD_FILE = 'IMOS_ANMN-NRS_BAD_VELOCITY_FILE.nc'
INPUT_FILES = [
    'IMOS_ANMN-NRS_AETVZ_20181213T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1812-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20181215T100000Z_C-20200430T000000Z.nc',
    'IMOS_ANMN-NRS_AETVZ_20180816T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1808-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20180822T053000Z_C-20200623T000000Z.nc',
    'IMOS_ANMN-NRS_AETVZ_20191016T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1910-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20191018T100000Z_C-20200430T000000Z.nc',
    BAD_FILE
]

OBS_VARS = {'TIME', 'DEPTH', 'DEPTH_quality_control', 'UCUR', 'UCUR_quality_control',
            'VCUR', 'VCUR_quality_control', 'WCUR', 'WCUR_quality_control', 'instrument_index', 'CELL_INDEX'}
INST_VARS = {'LATITUDE', 'LONGITUDE', 'NOMINAL_DEPTH', 'SECONDS_TO_MIDDLE'}
STR_VARS = {'source_file', 'instrument_id'}


class TestVelocityAggregatedTimeseries(BaseTestCase):
    EXPECTED_OUTPUT_FILE = os.path.join(
        TEST_ROOT, 'IMOS_ANMN-NRS_VZ_20180816_NRSROT_FV01_velocity-aggregated-timeseries_END-20191018_C-20200623.nc'
    )

    def test_velocity_aggregated(self):
        output_file, bad_files = velocity_aggregated(INPUT_FILES, 'NRSROT', input_dir=TEST_ROOT, output_dir='/tmp')

        self.assertEqual(1, len(bad_files))
        for file, errors in bad_files.items():
            self.assertEqual(BAD_FILE, file)

        dataset = Dataset(output_file)

        # check dimensions and variables
        self.assertSetEqual(set(dataset.dimensions), {'OBSERVATION', 'INSTRUMENT', 'strlen'})
        self.assertSetEqual(set(dataset.variables.keys()), OBS_VARS | INST_VARS | STR_VARS)

        for var in OBS_VARS:
            self.assertEqual(dataset.variables[var].dimensions, ('OBSERVATION',))
        for var in INST_VARS:
            self.assertEqual(dataset.variables[var].dimensions, ('INSTRUMENT',))
        for var in STR_VARS:
            self.assertEqual(dataset.variables[var].dimensions, ('INSTRUMENT', 'strlen'))

        for f in chartostring(dataset['source_file'][:]):
            self.assertIn(f, INPUT_FILES)

        # check attributes
        self.assertEqual(__version__, dataset.generating_code_version)
        self.assertIn(__version__, dataset.lineage)

        self.compare_global_attributes(dataset)

        self.check_nan_values(dataset)

        self.compare_variables(dataset)

    def test_all_rejected(self):
        self.assertRaises(NoInputFilesError, velocity_aggregated, [BAD_FILE], 'NRSROT',
                          input_dir=TEST_ROOT, output_dir='/tmp')


if __name__ == '__main__':
    unittest.main()
