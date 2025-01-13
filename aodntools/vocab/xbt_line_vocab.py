"""
Parse XBT line vocabulary from vocabs.ands.org.au
"""

import os
import urllib.request, urllib.error, urllib.parse
import xml.etree.ElementTree as ET
try:
    from functools import lru_cache
except ImportError:  # pragma: no cover
    from functools32 import lru_cache

DEFAULT_XBT_LINE_VOCAB_URL = 'http://content.aodn.org.au/Vocabularies/XBT-line/aodn_aodn-xbt-line-vocabulary.rdf'

class XbtLineVocabHelper(object):
    def __init__(self, xbt_line_vocab_url):
        self.xbt_line_vocab_url = xbt_line_vocab_url

    @lru_cache(maxsize=32)
    def xbt_line_info(self):
        """
        retrieves a dictionary of xbt line codes with their IMOS code equivalent if available
        """
        response = urllib.request.urlopen(self.xbt_line_vocab_url)
        html = response.read()
        root = ET.fromstring(html)

        xbt_dict = {}

        for item in root:
            if 'Description' in item.tag:
                xbt_line_code = None  # xbt_line_code
                xbt_line_pref_label = None  # xbt_line_code IMOS preferred value
                xbt_line_alt_label = None  # xbt line description, sometimes multiple altlabels such as AX01,
                # AX10, but these are not used by IMOS/AODN, so no need to make the code more complicated

                for val in item:
                    platform_element_sublabels = val.tag
                    if platform_element_sublabels is not None:
                        if 'prefLabel' in platform_element_sublabels:
                            xbt_line_pref_label = val.text
                        if 'code' in platform_element_sublabels:
                            xbt_line_code = val.text
                        if 'altLabel' in platform_element_sublabels:  #
                            xbt_line_alt_label = val.text

                if xbt_line_code is None and xbt_line_pref_label is not None:
                    xbt_dict[xbt_line_pref_label] = ({
                        'xbt_pref_label': xbt_line_pref_label,
                        'xbt_line_description': xbt_line_alt_label
                    })
                elif xbt_line_code is not None:
                    xbt_dict[xbt_line_code] = ({
                        'xbt_pref_label': xbt_line_pref_label,
                        'xbt_line_description': xbt_line_alt_label
                    })

        response.close()
        return xbt_dict