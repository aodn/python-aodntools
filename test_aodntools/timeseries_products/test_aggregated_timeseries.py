import unittest

from test_aodntools.base_test import BaseTestCase


class MyTestCase(BaseTestCase):
    def test_something(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
