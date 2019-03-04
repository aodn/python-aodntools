import unittest

from datetime import datetime, timedelta

from ncwriter.template import DatasetTemplate
from ncwriter.imos_template import ImosTemplate
from test_ncwriter.test_template import BaseTestCase


class TestImosTemplate(BaseTestCase):

    def setUp(self):
        super(TestImosTemplate, self).setUp()
        self.template = ImosTemplate()

    def test_init(self):
        self.assertIsInstance(self.template, DatasetTemplate)

    def test_fixed_global_attributes(self):
        global_attr = {
            "Conventions": "CF-1.6,IMOS-1.4",
            "project": "Integrated Marine Observing System (IMOS)",
            "naming_authority": "IMOS",
            "data_centre": "Australian Ocean Data Network (AODN)",
            "data_centre_email": "info@aodn.org.au"
        }
        for name, value in global_attr.items():
            self.assertEqual(value, self.template.global_attributes[name])

    def test_date_created(self):
        now = datetime.now()
        self.assertIsInstance(self.template.date_created, datetime)
        self.assertTrue(self.template.date_created - now < timedelta(seconds=1))

    def test_add_date_created_attribute(self):
        self.template.add_date_created_attribute()
        parsed_date_created = datetime.strptime(self.template.global_attributes['date_created'], '%Y-%m-%dT%H:%M:%SZ')
        self.assertTrue(self.template.date_created - parsed_date_created < timedelta(seconds=1))


if __name__ == '__main__':
    unittest.main()
