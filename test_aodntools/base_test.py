import shutil
import tempfile
import unittest

import numpy as np
from netCDF4 import Dataset


class BaseTestCase(unittest.TestCase):
    EXPECTED_OUTPUT_FILE = None

    @property
    def temp_dir(self):
        if not hasattr(self, '_temp_dir'):
            self._temp_dir = tempfile.mkdtemp(prefix=self.__class__.__name__)
        return self._temp_dir

    @property
    def temp_nc_file(self):
        if not hasattr(self, '_temp_nc_file'):
            with tempfile.NamedTemporaryFile(suffix='.nc', prefix=self.__class__.__name__, dir=self.temp_dir) as f:
                pass
            self._temp_nc_file = f.name
        return self._temp_nc_file

    def tearDown(self):
        if hasattr(self, '_temp_dir'):
            shutil.rmtree(self._temp_dir)

    def compare_global_attributes(self, dataset,
                                  attrs = ('geospatial_lat_max', 'geospatial_lat_min',
                                           'geospatial_lon_max', 'geospatial_lon_min',
                                           'geospatial_vertical_max', 'geospatial_vertical_min',
                                           'time_coverage_start', 'time_coverage_end'
                                           )
                                  ):
        "Compare global attributes of the given dataset with those in self.EXPECTED_OUTPUT_FILE"

        not_matching = []
        with Dataset(self.EXPECTED_OUTPUT_FILE) as expected:
            for attr in attrs:
                if dataset.getncattr(attr) != expected.getncattr(attr):
                    not_matching.append((attr,
                                         "expected: {exp}; found: {found}".format(exp=dataset.getncattr(attr),
                                                                                  found=dataset.getncattr(attr))
                                         ))

        self.assertEqual([], not_matching)

    def check_nan_values(self, dataset):
        "check that there are no NaN values in any variable (they should be fill values instead)"
        nan_vars = [(name, "contains NaN values")
                    for name, var in dataset.variables.items()
                    if var.dtype in (np.dtype('float32'), np.dtype('float64')) and any(np.isnan(var[:]))
                    ]
        self.assertEqual([], nan_vars)

    def compare_variables(self, dataset, skip_vars=('source_file', 'instrument_id')):
        """Compare dimensions and values of all variables in dataset with those in self.EXPECTED_OUTPUT_FILE,
        except for variables listed in skip_vars.
        """
        differences = []
        with Dataset(self.EXPECTED_OUTPUT_FILE) as expected:
            for var in set(expected.variables.keys()) - set(skip_vars):
                if not dataset[var].dimensions == expected[var].dimensions:
                    differences.append((var, "dimensions differ"))
                if not dataset[var].shape == expected[var].shape:
                    differences.append((var, "shapes differ"))

                # compare the raw data arrays (not the masked_array)
                if not all(np.isclose(dataset[var][:].data, expected[var][:].data)):
                    differences.append((var, "variable values differ"))

        self.assertEqual([], differences)
