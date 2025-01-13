import os
from unittest.mock import patch

from aodncore.pipeline import HandlerBase
from aodncore.testlib import BaseTestCase, HandlerTestCase
from aodndata.vocab.platform_code_vocab import (PlatformVocabHelper, platform_type_uris_by_category,
                                                platform_altlabels_per_preflabel)

TEST_ROOT = os.path.join(os.path.dirname(__file__))

TEST_PLATFORM_CAT_VOCAB_URL = '%s%s' % ('file://',
                                        os.path.join(TEST_ROOT, 'aodn_aodn-platform-category-vocabulary.rdf'))

TEST_PLATFORM_VOCAB_URL = '%s%s' % ('file://',
                                    os.path.join(TEST_ROOT, 'aodn_aodn-platform-vocabulary.rdf'))


class TestVocabHandler(HandlerBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.platform_vocab_helper = PlatformVocabHelper(TEST_PLATFORM_VOCAB_URL, TEST_PLATFORM_CAT_VOCAB_URL)


class TestDummyHandler(HandlerTestCase):
    def setUp(self):
        self.handler_class = TestVocabHandler
        super().setUp()

    def test_platform_vocab_helper(self):
        handler = self.handler_class(self.temp_nc_file)
        self.assertIsInstance(handler.platform_vocab_helper.platform_type_uris_by_category(), dict)
        self.assertEqual(handler.platform_vocab_helper.platform_altlabels_per_preflabel('Vessel')['9VUU'], 'Anro Asia')


class TestPlatformCodeVocab(BaseTestCase):
    def setUp(self):
        self.platform_vocab_helper = PlatformVocabHelper(TEST_PLATFORM_VOCAB_URL, TEST_PLATFORM_CAT_VOCAB_URL)

    def test_platform_type(self):
        res = self.platform_vocab_helper.platform_type_uris_by_category()
        self.assertEqual(res['Float'][0], 'http://vocab.nerc.ac.uk/collection/L06/current/42')

    def test_altlabels(self):
        res = self.platform_vocab_helper.platform_altlabels_per_preflabel(category_name='Vessel')
        self.assertEqual(res['9VUU'], 'Anro Asia')

    def test_platform_type_deprecated(self):
        res = platform_type_uris_by_category(platform_cat_vocab_url=TEST_PLATFORM_CAT_VOCAB_URL)
        self.assertEqual(res['Float'][0], 'http://vocab.nerc.ac.uk/collection/L06/current/42')

    @patch("aodndata.vocab.platform_code_vocab.DEFAULT_PLATFORM_CAT_VOCAB_URL", TEST_PLATFORM_CAT_VOCAB_URL)
    def test_altlabels_deprecated(self):
        res = platform_altlabels_per_preflabel(category_name='Vessel', platform_vocab_url=TEST_PLATFORM_VOCAB_URL)
        self.assertEqual(res['9VUU'], 'Anro Asia')
