# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.sdk import PluginManager

from q2_meta_vsearch.plugin_setup import plugin as vsearch_plugin


class PluginSetupTests(unittest.TestCase):

    def test_plugin_setup(self):
        self.assertEqual(vsearch_plugin.name, 'meta-vsearch')

    def test_coexists_with_vsearch(self):
        plugins = PluginManager().plugins
        self.assertIn('vsearch', plugins)
        self.assertIn('meta-vsearch', plugins)
        upstream = plugins['vsearch']
        self.assertTrue(set(upstream.actions).issubset(vsearch_plugin.actions))
        self.assertFalse(set(upstream.formats) & set(vsearch_plugin.formats))
        self.assertFalse(set(upstream.type_fragments) &
                         set(vsearch_plugin.type_fragments))
