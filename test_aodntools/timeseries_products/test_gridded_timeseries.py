#!/usr/bin/env python3

import os
import unittest

from netCDF4 import Dataset

from test_aodntools.base_test import BaseTestCase
from aodntools import __version__
from aodntools.timeseries_products.gridded_timeseries import grid_variable


TEST_ROOT = os.path.dirname(__file__)
INPUT_FILE = 'IMOS_ANMN-NRS_STZ_20181213_NRSROT_FV02_hourly-timeseries_END-20190523_C-20220428.nc'


class TestGriddedTimeseries(BaseTestCase):
    EXPECTED_OUTPUT_FILE = os.path.join(
        TEST_ROOT, 'IMOS_ANMN-NRS_TZ_20181213_NRSROT_FV02_TEMP-gridded-timeseries_END-20190523_C-20230110.nc'
    )

    def test_grid_variable(self):
        output_file = grid_variable(INPUT_FILE, 'TEMP', input_dir=TEST_ROOT, output_dir='/tmp')

        self.assertRegex(output_file,
                         r'IMOS_ANMN-NRS_TZ_20181213_NRSROT_FV02_TEMP-gridded-timeseries_END-20190523_C-\d{8}\.nc'
                         )

        dataset = Dataset(output_file)
        self.assertSetEqual(set(dataset.dimensions), {'TIME', 'DEPTH'})
        self.assertSetEqual(set(dataset.variables.keys()),
                            {'TIME', 'DEPTH', 'LATITUDE', 'LONGITUDE', 'TEMP', 'TEMP_count'})

        # check metadata
        self.assertEqual(__version__, dataset.generating_code_version)
        self.assertIn(__version__, dataset.lineage)
        self.assertIn('gridded_timeseries.py', dataset.lineage)
        self.assertIn(INPUT_FILE, dataset.source_file)

        self.compare_global_attributes(dataset)

        self.check_nan_values(dataset)

        self.compare_variables(dataset)


if __name__ == '__main__':
    unittest.main()
