import os

from aodncore.testlib import BaseTestCase
from aodndata.vocab.xbt_line_vocab import XbtLineVocabHelper

TEST_ROOT = os.path.join(os.path.dirname(__file__))

TEST_XBT_LINE_VOCAB_URL = '%s%s' % ('file://',
                                    os.path.join(TEST_ROOT, 'aodn_aodn-xbt-line-vocabulary.rdf'))


class TestPlatformCodeVocab(BaseTestCase):
    def setUp(self):

        self.xbt_line_vocab_helper = XbtLineVocabHelper(TEST_XBT_LINE_VOCAB_URL)

    def test_platform_type(self):
        res = self.xbt_line_vocab_helper.xbt_line_info()

        self.assertTrue('AX01' in res.keys())
        self.assertEqual('Greenland - Iceland - Scotland', res['AX01']['xbt_line_description'])
