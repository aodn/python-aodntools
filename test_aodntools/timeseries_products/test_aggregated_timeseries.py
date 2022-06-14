#!/usr/bin/env python3

import os
import unittest

from netCDF4 import Dataset, chartostring

from aodntools import __version__
from aodntools.timeseries_products.aggregated_timeseries import main_aggregator
from aodntools.timeseries_products.common import NoInputFilesError
from test_aodntools.base_test import BaseTestCase

TEST_ROOT = os.path.dirname(__file__)
BAD_FILE = 'IMOS_ANMN-NRS_TZ_20181213T080000Z_NRSROT_FV00_NRSROT-1812-SBE39-43_END-20181214T004000Z_C-20190827T000000Z.nc'
INPUT_FILES = [
    'IMOS_ANMN-NRS_TZ_20181213T080000Z_NRSROT_FV01_NRSROT-1812-SBE39-23_END-20190306T160000Z_C-20190827T000000Z.nc',
    'IMOS_ANMN-NRS_TZ_20190313T144000Z_NRSROT_FV01_NRSROT-1903-SBE39-27_END-20190524T010000Z_C-20190827T000000Z.nc',
    'IMOS_ANMN-NRS_BCKOSTUZ_20181213T080038Z_NRSROT_FV01_NRSROT-1812-WQM-55_END-20181215T013118Z_C-20190828T000000Z.nc',
    BAD_FILE
]


class TestAggregatedTimeseries(BaseTestCase):
    EXPECTED_OUTPUT_FILE = os.path.join(
        TEST_ROOT, 'IMOS_ANMN-NRS_TZ_20181213_NRSROT_FV01_TEMP-aggregated-timeseries_END-20190523_C-20220607.nc'
    )

    def test_main_aggregator(self):
        output_file, bad_files = main_aggregator(INPUT_FILES, 'TEMP', 'NRSROT', input_dir=TEST_ROOT,
                                                 output_dir='/tmp')

        self.assertEqual(1, len(bad_files))
        for file, errors in bad_files.items():
            self.assertEqual(BAD_FILE, file)
            self.assertSetEqual(set(errors), {'no NOMINAL_DEPTH',
                                              'Wrong file version: Level 0 - Raw Data',
                                              'no time_deployment_start attribute',
                                              'no time_deployment_end attribute'
                                              }
                                )

        dataset = Dataset(output_file)

        # check dimensions and variables
        self.assertSetEqual(set(dataset.dimensions), {'OBSERVATION', 'INSTRUMENT', 'strlen'})
        self.assertSetEqual(set(dataset.variables.keys()),
                            {'TIME', 'LATITUDE', 'LONGITUDE', 'NOMINAL_DEPTH', 'DEPTH', 'DEPTH_quality_control',
                             'PRES', 'PRES_quality_control', 'PRES_REL', 'PRES_REL_quality_control',
                             'TEMP', 'TEMP_quality_control', 'instrument_index', 'instrument_id', 'source_file'}
                            )

        obs_vars = {'TIME', 'DEPTH', 'DEPTH_quality_control', 'PRES', 'PRES_quality_control',
                    'PRES_REL', 'PRES_REL_quality_control', 'TEMP', 'TEMP_quality_control', 'instrument_index'}
        for var in obs_vars:
            self.assertEqual(dataset.variables[var].dimensions, ('OBSERVATION',))

        inst_vars = {'LATITUDE', 'LONGITUDE', 'NOMINAL_DEPTH'}
        for var in inst_vars:
            self.assertEqual(dataset.variables[var].dimensions, ('INSTRUMENT',))

        string_vars = {'source_file', 'instrument_id'}
        for var in string_vars:
            self.assertEqual(dataset.variables[var].dimensions, ('INSTRUMENT', 'strlen'))

        for f in chartostring(dataset['source_file'][:]):
            self.assertIn(f, INPUT_FILES)

        # check attributes
        self.assertEqual(__version__, dataset.generating_code_version)
        self.assertIn(__version__, dataset.lineage)
        self.assertIn(BAD_FILE, dataset.rejected_files)

        self.compare_global_attributes(dataset)

        self.check_nan_values(dataset)

        self.compare_variables(dataset)

    def test_source_file_attributes(self):
        output_file, bad_files = main_aggregator(INPUT_FILES, 'PSAL', 'NRSROT', input_dir=TEST_ROOT,
                                                 output_dir='/tmp', download_url_prefix='http://test.download.url',
                                                 opendap_url_prefix='http://test.opendap.url'
                                                 )
        dataset = Dataset(output_file)
        self.assertEqual(dataset['source_file'].download_url_prefix, 'http://test.download.url')
        self.assertEqual(dataset['source_file'].opendap_url_prefix, 'http://test.opendap.url')

    def test_all_rejected(self):
        self.assertRaises(NoInputFilesError, main_aggregator, [BAD_FILE], 'TEMP', 'NRSROT',
                          input_dir=TEST_ROOT, output_dir='/tmp')


if __name__ == '__main__':
    unittest.main()
