#!/usr/bin/env python3

import os
import unittest

from netCDF4 import Dataset, chartostring

from test_aodntools.base_test import BaseTestCase
from aodntools import __version__
from aodntools.timeseries_products.common import NoInputFilesError
from aodntools.timeseries_products.hourly_timeseries import hourly_aggregator


TEST_ROOT = os.path.dirname(__file__)
BAD_FILE = 'IMOS_ANMN-NRS_TZ_20181213T080000Z_NRSROT_FV00_NRSROT-1812-SBE39-43_END-20181214T004000Z_C-20190827T000000Z.nc'
INPUT_FILES = [
    'IMOS_ANMN-NRS_BCKOSTUZ_20181213T080038Z_NRSROT_FV01_NRSROT-1812-WQM-55_END-20181215T013118Z_C-20190828T000000Z.nc',
    'IMOS_ANMN-NRS_TZ_20181213T080000Z_NRSROT_FV01_NRSROT-1812-SBE39-23_END-20190306T160000Z_C-20190827T000000Z.nc',
    'IMOS_ANMN-NRS_TZ_20190313T144000Z_NRSROT_FV01_NRSROT-1903-SBE39-27_END-20190524T010000Z_C-20190827T000000Z.nc',
    BAD_FILE
]
INPUT_PATHS = [os.path.join(TEST_ROOT, f) for f in INPUT_FILES]

INST_VARIABLES = {'instrument_id', 'source_file', 'LONGITUDE', 'LATITUDE', 'NOMINAL_DEPTH'}
OBS_VARIABLES = {'instrument_index', 'TIME'}
measured_variables = {'DEPTH', 'CPHL', 'CHLF', 'CHLU',
                      'DOX', 'DOX1', 'DOX1_2', 'DOX1_3', 'DOX2', 'DOX2_1', 'DOXS', 'DOXY',
                      'PRES', 'PRES_REL', 'PSAL', 'TEMP', 'TURB', 'PAR'
                      }
function_stats = ['_min', '_max', '_std', '_count']
for v in measured_variables:
    OBS_VARIABLES.add(v)
    for s in function_stats:
        OBS_VARIABLES.add(v + s)

NO_INWATER_DATA_FILE = 'IMOS_ANMN-NSW_TZ_PH100_NO_INWATER_DATA.nc'
PH100_FILES = [
    'IMOS_ANMN-NSW_TZ_20200703T001500Z_PH100_FV01_PH100-2007-Aqualogger-520T-96_END-20200907T233000Z_C-20210112T044909Z.nc',
    'IMOS_ANMN-NSW_TZ_PH100_ALL_FLAGGED_BAD.nc',
    NO_INWATER_DATA_FILE
]

SYD100_FILES = [
    'IMOS_ANMN-NSW_TZ_SYD100_BAD_TIMESTAMPS.nc',
]


class TestHourlyTimeseries(BaseTestCase):
    EXPECTED_OUTPUT_FILE = os.path.join(
        TEST_ROOT, 'IMOS_ANMN-NRS_STZ_20181213_NRSROT_FV02_hourly-timeseries_END-20190523_C-20220428.nc'
    )

    def test_hourly_aggregator(self):
        output_file, bad_files = hourly_aggregator(files_to_aggregate=INPUT_PATHS,
                                                   site_code='NRSROT',
                                                   qcflags=(1, 2),
                                                   output_dir='/tmp'
                                                   )
        self.assertRegex(output_file,
                         r'IMOS_ANMN-NRS_STZ_20181213_NRSROT_FV02_hourly-timeseries_END-20190523_C-\d{8}\.nc'
                         )

        self.assertEqual(1, len(bad_files))
        for path, errors in bad_files.items():
            self.assertEqual(os.path.join(TEST_ROOT, BAD_FILE), path)
            self.assertSetEqual(set(errors), {'no NOMINAL_DEPTH',
                                              'Wrong file version: Level 0 - Raw Data',
                                              'no time_deployment_start attribute',
                                              'no time_deployment_end attribute'
                                              }
                                )

        dataset = Dataset(output_file)
        self.assertSetEqual(set(dataset.dimensions), {'OBSERVATION', 'INSTRUMENT', 'string256'})
        self.assertTrue(set(dataset.variables.keys()).issubset(OBS_VARIABLES | INST_VARIABLES))

        inst_variables = {n for n, v in dataset.variables.items() if v.dimensions[0] == 'INSTRUMENT'}
        self.assertSetEqual(inst_variables, INST_VARIABLES)
        obs_variables = {n for n, v in dataset.variables.items() if v.dimensions == ('OBSERVATION',)}
        self.assertTrue(obs_variables.issubset(OBS_VARIABLES))

        for f in chartostring(dataset['source_file'][:]):
            self.assertIn(f, INPUT_PATHS)

        # check metadata
        self.assertEqual(__version__, dataset.generating_code_version)
        self.assertIn(__version__, dataset.lineage)
        self.assertIn('hourly_timeseries.py', dataset.lineage)
        self.assertIn(BAD_FILE, dataset.rejected_files)

        self.compare_global_attributes(dataset)

        self.check_nan_values(dataset)

        self.compare_variables(dataset)

    def test_hourly_aggregator_with_nonqc(self):
        output_file, bad_files = hourly_aggregator(files_to_aggregate=INPUT_FILES,
                                                   site_code='NRSROT',
                                                   qcflags=(0, 1, 2),
                                                   input_dir=TEST_ROOT,
                                                   output_dir='/tmp',
                                                   download_url_prefix='http://test.download.url',
                                                   opendap_url_prefix='http://test.opendap.url'
                                                   )
        self.assertRegex(output_file,
                         r'IMOS_ANMN-NRS_BOSTUZ_20181213_NRSROT_FV02_hourly-timeseries-including-non-QC'
                         r'_END-20190523_C-\d{8}\.nc'
                         )

        dataset = Dataset(output_file)
        self.assertEqual(dataset['source_file'].download_url_prefix, 'http://test.download.url')
        self.assertEqual(dataset['source_file'].opendap_url_prefix, 'http://test.opendap.url')
        for f in chartostring(dataset['source_file'][:]):
            self.assertIn(f, INPUT_FILES)

    def test_with_adcp(self):
        # Replace the BAD_FILE with an ADCP file - aggregation should work (only takes TEMP from the ADCP)
        input_files = INPUT_FILES[:2] + \
                      ['IMOS_ANMN-NRS_AETVZ_20180816T080000Z_NRSROT-ADCP_FV01_NRSROT-ADCP-1808-Sentinel-or-Monitor-Workhorse-ADCP-44_END-20180822T053000Z_C-20200623T000000Z.nc']
        output_file, bad_files = hourly_aggregator(files_to_aggregate=input_files,
                                                   site_code='NRSROT',
                                                   qcflags=(1, 2),
                                                   input_dir=TEST_ROOT,
                                                   output_dir='/tmp'
                                                   )

        self.assertEqual(0, len(bad_files))

    def test_all_rejected(self):
        self.assertRaises(NoInputFilesError, hourly_aggregator, [BAD_FILE], 'NRSROT', (1, 2), input_dir=TEST_ROOT)

    def test_some_files_without_good_data(self):
        output_file, bad_files = hourly_aggregator(files_to_aggregate=PH100_FILES,
                                                   site_code='PH100',
                                                   qcflags=(1, 2),
                                                   input_dir=TEST_ROOT,
                                                   output_dir='/tmp'
                                                   )
        self.assertEqual(1, len(bad_files))
        for path, errors in bad_files.items():
            self.assertEqual(NO_INWATER_DATA_FILE, path)
            self.assertIn('no in-water data', errors)
        with Dataset(output_file) as dataset:
            self.check_nan_values(dataset)

    def test_bad_timestamps(self):
        output_file, bad_files = hourly_aggregator(files_to_aggregate=SYD100_FILES,
                                                   site_code='SYD100',
                                                   qcflags=(1, 2),
                                                   input_dir=TEST_ROOT,
                                                   output_dir='/tmp'
                                                   )
        self.assertEqual(0, len(bad_files))


if __name__ == '__main__':
    unittest.main()
