import shutil
import tempfile
import unittest


class BaseTestCase(unittest.TestCase):

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
