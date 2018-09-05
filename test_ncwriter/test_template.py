import unittest

from ncwriter import DatasetTemplate


class TestDatasetTemplate(unittest.TestCase):
    def test_init(self):
        template = DatasetTemplate()

# Functionality to test
#
# TODO: create template from dicts
# e.g. template = DatasetTemplate(dims, vars, globalattrs)
# TODO: create template from json file
# e.g. DatasetTemplate.from_json(path='template.json')
# TODO: create template from other formats (later...)

# TODO: add global attributes
# e.g. template.title = 'Test dataset'

# TODO: add dimensions
# e.g. template.dimensions['TIME'] = 100

# TODO: add variables
# TODO: add variable attributes
# e.g. template.variables['PRES']['units'] = 'dbar'

# TODO: add data from numpy arrays
# e.g. template.add_data(TIME=time_values, TEMP=temp_values, PRES=pres_values)
# TODO: add data from Pandas dataframe (later...)
# e.g. template.add_data(dataframe)

# TODO: create netCDF file
# e.g. template.create(filename)   # user-specified file name
# TODO: create netCDF file with auto-generated file name according to IMOS conventions
# e.g. template.create()


if __name__ == '__main__':
    unittest.main()
