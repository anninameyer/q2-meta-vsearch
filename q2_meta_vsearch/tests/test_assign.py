# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

from fractions import Fraction
import json
from pathlib import Path
import platform
from unittest.mock import patch

import pandas as pd
from qiime2 import Artifact, Metadata
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (
    BLAST6Format, SequenceCharacteristicsDirectoryFormat, TSVTaxonomyFormat)

from q2_meta_vsearch._assign import read_reference_taxonomy
from q2_meta_vsearch._assignment import (
    assign_reference_alignments, SUMMARY_COLUMNS)
from q2_meta_vsearch._format import MetaVsearchReferenceAlignmentsDirFmt
from q2_meta_vsearch._taxonomy import parse_lineage


class AssignmentTests(TestPluginBase):
    package = 'q2_meta_vsearch.tests'

    def setUp(self):
        super().setUp()
        # Newer framework provenance may call the OS `file` utility here.
        # Cache platform information before guarding assignment execution.
        platform.platform()
        self.external = patch('subprocess.run', side_effect=AssertionError(
            'Assignment must not execute external commands')).start()
        self.addCleanup(patch.stopall)

    def _path(self, name):
        return Path(self.get_data_path('assignment/' + name))

    def _assign(self, taxonomy='two-species.tsv', evidence='evidence',
                **kwargs):
        return assign_reference_alignments(
            self._path(evidence),
            read_reference_taxonomy(self._path(taxonomy)),
            **kwargs)

    def _row(self, result, query='00123'):
        return next(row for row in result.summary if row['id'] == query)

    def _action(self, taxonomy='two-species.tsv', evidence='evidence',
                **kwargs):
        alignments = Artifact.import_data(
            'FeatureData[ReferenceAlignments]', str(self._path(evidence)),
            view_type=MetaVsearchReferenceAlignmentsDirFmt)
        taxa = Artifact.import_data('FeatureData[Taxonomy]',
                                    str(self._path(taxonomy)),
                                    view_type=TSVTaxonomyFormat)
        return self.plugin.methods['assign_closed_reference'](
            alignments=alignments, reference_taxonomy=taxa, **kwargs)

    def test_unique_reference(self):
        result = self._assign(evidence='evidence-unique')
        row = self._row(result, 'unique')
        self.assertEqual(row['n_best_reference_hits'], 1)
        self.assertEqual(row['reference_without_weights'], 'ref-a')
        self.assertEqual(row['reference_assignment_status'],
                         'unique_reference')
        self.assertEqual(row['taxonomic_resolution_rank'], 'Species')
        self.assertEqual(result.mapping, {'unique': 'ref-a'})

    def test_same_species_reference_ambiguity(self):
        row = self._row(self._assign('same-species.tsv'))
        self.assertEqual(row['n_best_reference_hits'], 3)
        self.assertEqual(row['n_best_alignments'], 3)
        self.assertEqual(row['n_best_taxa'], 1)
        self.assertEqual(row['taxonomic_assignment_status'], 'resolved')
        self.assertEqual(row['taxonomic_resolution_rank'], 'Species')
        self.assertEqual(row['taxonomy_ambiguity_label'],
                         'Lactiplantibacillus plantarum')
        self.assertEqual(row['reference_assignment_status'],
                         'deterministic_reference')

    def test_two_species_lca_and_label(self):
        row = self._row(self._assign())
        self.assertEqual(row['n_best_reference_hits'], 3)
        self.assertEqual(row['n_best_taxa'], 2)
        self.assertEqual(row['taxonomic_resolution_rank'], 'Genus')
        self.assertEqual(row['taxonomic_assignment_status'], 'ambiguous')
        self.assertTrue(row['taxonomy_lca'].endswith(
            'g__Lactiplantibacillus'))
        self.assertNotIn('s__', row['taxonomy_without_weights'])
        self.assertEqual(row['taxonomy_ambiguity_label'],
                         'Lactiplantibacillus paraplantarum|plantarum')

    def test_compact_label_limit(self):
        row = self._row(self._assign('three-species.tsv'))
        self.assertEqual(row['n_best_taxa'], 3)
        self.assertEqual(row['taxonomy_ambiguity_label'],
                         'Lactiplantibacillus ssg.')
        explicit = self._row(self._assign('three-species.tsv',
                                          max_explicit_species=3))
        self.assertNotIn('ssg.', explicit['taxonomy_ambiguity_label'])
        self.assertEqual(explicit['taxonomy_ambiguity_label'].count('|'), 2)
        self.assertEqual(row['taxonomy_lca'], explicit['taxonomy_lca'])
        self.assertNotIn('ssg.', row['taxonomy_without_weights'])

    def test_across_genera(self):
        row = self._row(self._assign('across-genera.tsv'))
        self.assertEqual(row['taxonomic_resolution_rank'], 'Family')
        self.assertEqual(row['taxonomy_ambiguity_label'], 'Lactobacillaceae')
        self.assertNotIn('g__', row['taxonomy_lca'])

    def test_across_higher_ranks(self):
        for taxonomy, rank, lca in [
                ('across-phyla.tsv', 'Domain', 'd__Bacteria'),
                ('across-domains.tsv', 'Unassigned', 'Unassigned')]:
            with self.subTest(taxonomy=taxonomy):
                row = self._row(self._assign(taxonomy))
                self.assertEqual(row['taxonomic_resolution_rank'], rank)
                self.assertEqual(row['taxonomy_lca'], lca)

    def test_explicit_kingdom_is_not_inferred_as_domain(self):
        row = self._row(self._assign('fungal-taxonomy.tsv'))
        self.assertEqual(row['taxonomic_resolution_rank'], 'Species')
        self.assertTrue(row['taxonomy_lca'].startswith('k__Fungi;'))
        kingdom = self._row(self._assign('kingdom-taxonomy.tsv'))
        self.assertEqual(kingdom['taxonomic_resolution_rank'], 'Kingdom')
        self.assertEqual(kingdom['taxonomy_lca'], 'k__Fungi')

    def test_missing_species(self):
        row = self._row(self._assign('missing-species.tsv'))
        self.assertEqual(row['taxonomic_assignment_status'], 'incomplete')
        self.assertEqual(row['taxonomic_resolution_rank'], 'Genus')
        self.assertEqual(row['n_best_taxa'], 1)
        self.assertEqual(row['taxonomy_ambiguity_label'],
                         'Lactiplantibacillus')
        self.assertNotIn('s__', row['taxonomy_lca'])

    def test_incomplete_depths(self):
        row = self._row(self._assign('incomplete.tsv'))
        self.assertEqual(row['taxonomic_resolution_rank'], 'Family')
        self.assertEqual(row['taxonomic_assignment_status'], 'ambiguous')

    def test_missing_reference_annotation(self):
        row = self._row(self._assign('missing-taxonomy.tsv'))
        self.assertEqual(row['n_missing_taxonomy_reference_hits'], 1)
        self.assertEqual(row['taxonomic_assignment_status'],
                         'partial_missing_taxonomy')
        self.assertEqual(row['taxonomy_lca'], 'Unassigned')
        self.assertEqual(row['n_best_taxa'], 1)
        self.assertEqual(row['reference_without_weights'], 'ref-a')
        excluded = self._row(self._assign('missing-taxonomy.tsv',
                                          unresolved_policy='exclude'))
        self.assertIsNone(excluded['reference_without_weights'])

    def test_all_reference_ids_absent_from_taxonomy(self):
        row = self._row(self._assign('all-missing.tsv'))
        self.assertEqual(row['n_missing_taxonomy_reference_hits'], 3)
        self.assertEqual(row['taxonomic_assignment_status'],
                         'missing_taxonomy')
        self.assertEqual(row['n_best_taxa'], 0)
        self.assertEqual(row['taxonomy_lca'], 'Unassigned')

    def test_nonstandard_species_are_conservative(self):
        row = self._row(self._assign('nonstandard-species.tsv'))
        self.assertEqual(row['taxonomic_resolution_rank'], 'Genus')
        self.assertEqual(row['taxonomy_ambiguity_label'],
                         'Lactiplantibacillus')
        self.assertNotIn('ssg.', row['taxonomy_lca'])

    def test_reference_multiplicity_is_not_a_vote(self):
        first = self._row(self._assign())
        second = self._row(self._assign(evidence='evidence-less-duplicates'))
        for key in ['taxonomy_lca', 'taxonomy_without_weights', 'n_best_taxa',
                    'taxonomic_assignment_status', 'taxonomy_ambiguity_label']:
            self.assertEqual(first[key], second[key])
        self.assertEqual(first['n_best_reference_hits'], 3)
        self.assertEqual(second['n_best_reference_hits'], 2)

    def test_reference_order_independence_and_selection(self):
        first = self._assign()
        second = self._assign(evidence='evidence-reversed')
        self.assertEqual(first, second)
        self.assertEqual(self._row(first)['reference_without_weights'],
                         'ref-a')

    def test_exclude_only_unresolved_taxonomy(self):
        result = self._assign(unresolved_policy='exclude')
        row = self._row(result)
        self.assertEqual(row['matched'], 'true')
        self.assertIsNone(row['reference_without_weights'])
        self.assertEqual(row['statistics_reference'], 'ref-a')
        self.assertEqual(row['n_remaining_reference_hits'], 3)
        self.assertEqual(row['reference_assignment_status'],
                         'excluded_taxonomic_ambiguity')
        self.assertNotIn('00123', result.mapping)
        self.assertEqual(row['best_identity'], 1)
        same = self._row(self._assign('same-species.tsv',
                                      unresolved_policy='exclude'))
        self.assertEqual(same['reference_without_weights'], 'ref-a')

    def test_exact_best_and_near_best_are_separate(self):
        base = self._row(self._assign())
        diagnostic = self._row(self._assign(identity_tolerance=1))
        self.assertEqual(base['n_near_best_reference_hits'], 0)
        self.assertEqual(diagnostic['n_near_best_reference_hits'], 1)
        self.assertEqual(base['n_accepted_reference_hits'], 4)
        for key in ['reference_without_weights', 'taxonomy_without_weights',
                    'n_best_reference_hits', 'n_best_taxa']:
            self.assertEqual(base[key], diagnostic[key])
        self.assertNotIn('Archaea', diagnostic['taxonomy_lca'])

    def test_precision_and_tolerance_boundary(self):
        result = self._assign('precision-taxonomy.tsv', 'evidence-precision',
                              identity_tolerance=0.05)
        row = self._row(result, 'precision')
        self.assertEqual(row['reference_without_weights'], 'better')
        self.assertEqual(row['n_best_reference_hits'], 1)
        self.assertEqual(row['n_near_best_reference_hits'], 1)
        self.assertEqual(Fraction(row['best_identity_numerator'],
                                  row['best_identity_denominator']),
                         Fraction(2998, 3000))
        self.assertEqual(row['taxonomic_resolution_rank'], 'Species')

    def test_threshold_consistency_90_through_100(self):
        baseline = self._row(self._assign())
        for percent in range(90, 101):
            with self.subTest(percent=percent):
                result = self._assign(identity_threshold=percent / 100)
                row = self._row(result)
                self.assertEqual(row['n_best_reference_hits'], 3)
                self.assertEqual(row['reference_without_weights'],
                                 baseline['reference_without_weights'])
                self.assertEqual(row['taxonomy_lca'], baseline['taxonomy_lca'])
                boundary = self._row(result, 'boundary')
                self.assertEqual(boundary['matched'],
                                 'false' if percent == 100 else 'true')
        boundary = self._row(self._assign(identity_threshold=0.995),
                             'boundary')
        self.assertEqual(boundary['matched'], 'true')

    def test_unmatched_rows_and_original_identifiers(self):
        result = self._assign()
        self.assertEqual({row['id'] for row in result.summary},
                         {'00123', 'NA', '1e3', 'asv;size=7;odd',
                          'boundary', 'unmatched'})
        self.assertEqual(len(result.summary), 6)
        row = self._row(result, 'unmatched')
        self.assertEqual(row['matched'], 'false')
        for key in ['best_identity', 'best_alignment_length', 'best_gaps',
                    'reference_without_weights', 'statistics_reference']:
            self.assertIsNone(row[key])
        self.assertEqual(row['n_best_reference_hits'], 0)
        self.assertEqual(result.taxonomy['unmatched'], 'Unassigned')

    def test_statistics_use_named_reference_not_best_coverage(self):
        result = self._assign('statistics-taxonomy.tsv', 'evidence-statistics')
        row = self._row(result, 'statistics')
        self.assertEqual(row['statistics_reference'], 'ref-a')
        self.assertEqual(row['best_alignment_length'], 160)
        self.assertEqual(row['best_query_coverage'], 0.8)
        self.assertEqual(row['best_reference_coverage'], 1)
        excluded = self._row(self._assign(
            'statistics-taxonomy.tsv', 'evidence-statistics',
            unresolved_policy='exclude'), 'statistics')
        self.assertEqual(row['best_alignment_length'],
                         excluded['best_alignment_length'])

    def test_metadata_taxonomy_artifacts_and_roundtrip(self):
        outputs = self._action()
        self.assertEqual(str(outputs.assignment_metadata.type),
                         'FeatureData[SequenceCharacteristics]')
        outputs.assignment_metadata.validate(level='max')
        outputs.taxonomy.validate(level='max')
        md = outputs.assignment_metadata.view(Metadata)
        self.assertEqual(len(md.ids), 6)
        self.assertIn('00123', md.ids)
        self.assertIn('NA', md.ids)
        frame = outputs.assignment_metadata.view(pd.DataFrame)
        self.assertEqual(list(frame.columns), list(SUMMARY_COLUMNS))
        self.assertEqual(frame.loc['00123', 'reference_without_weights'],
                         'ref-a')
        self.assertEqual(md.columns['n_best_reference_hits'].type, 'numeric')
        output = Path(self.temp_dir.name) / 'metadata.tsv'
        md.save(str(output))
        restored = Metadata.load(str(output))
        self.assertEqual(set(restored.ids), set(md.ids))
        tax = outputs.taxonomy.view(pd.DataFrame)
        self.assertEqual(set(tax.index), set(md.ids))
        self.assertEqual(set(outputs.taxonomy.view(pd.Series).index),
                         set(md.ids))
        self.assertNotIn('s__', tax.loc['00123', 'Taxon'])
        self.assertNotIn('ssg.', tax.loc['00123', 'Taxon'])
        saved = outputs.assignment_metadata.save(
            str(Path(self.temp_dir.name) / 'assignments'))
        loaded = Artifact.load(saved)
        self.assertEqual(set(loaded.view(Metadata).ids), set(md.ids))
        exported = Path(self.temp_dir.name) / 'exported'
        outputs.assignment_metadata.export_data(str(exported))
        reimported = Artifact.import_data(
            'FeatureData[SequenceCharacteristics]', str(exported),
            view_type=SequenceCharacteristicsDirectoryFormat)
        self.assertEqual(set(reimported.view(Metadata).ids), set(md.ids))

    def test_numeric_only_ids_and_reference_ids(self):
        outputs = self._action('numeric-taxonomy.tsv', 'evidence-numeric')
        metadata = outputs.assignment_metadata.view(Metadata)
        self.assertEqual(set(metadata.ids), {'00123', '00045'})
        self.assertEqual(metadata.to_dataframe().loc[
            '00123', 'reference_without_weights'], '0001')
        self.assertEqual(set(outputs.taxonomy.view(pd.DataFrame).index),
                         {'00123', '00045'})

    def test_all_unmatched_metadata_types(self):
        outputs = self._action(evidence='evidence-unmatched')
        outputs.taxonomy.validate(level='max')
        metadata = outputs.assignment_metadata.view(Metadata)
        self.assertEqual(metadata.columns['best_identity'].type, 'numeric')
        self.assertTrue(metadata.to_dataframe()['best_identity'].isna().all())

    def test_missing_taxonomy_and_duplicate_ids(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            read_reference_taxonomy(self._path('duplicate-taxonomy.tsv'))

    def test_supported_lineages_and_missing_ranks(self):
        cases = json.loads(self._path('lineages.json').read_text())
        for case in cases:
            with self.subTest(lineage=case['input']):
                self.assertEqual(parse_lineage(case['input']).taxonomy,
                                 case['expected'])

    def test_invalid_parameters(self):
        for params in [{'identity_threshold': 0.89},
                       {'identity_threshold': 1.1},
                       {'identity_tolerance': -1}, {'identity_tolerance': 101},
                       {'unresolved_policy': 'random'},
                       {'max_explicit_species': 0}]:
            with self.subTest(parameters=params):
                with self.assertRaises(ValueError):
                    self._assign(**params)

    def test_no_vsearch_execution(self):
        self._assign()
        self._action()
        self.external.assert_not_called()

    def test_blast6_no_hit_records_are_not_assignment_input(self):
        blast6 = Artifact.import_data('FeatureData[BLAST6]',
                                      str(self._path('no-hit.blast6.tsv')),
                                      view_type=BLAST6Format)
        taxonomy = Artifact.import_data(
            'FeatureData[Taxonomy]', str(self._path('two-species.tsv')),
            view_type=TSVTaxonomyFormat)
        with self.assertRaises(TypeError):
            self.plugin.methods['assign_closed_reference'](
                alignments=blast6, reference_taxonomy=taxonomy)

    def test_no_weighted_placeholder_columns(self):
        result = self._assign()
        columns = set(self._row(result))
        self.assertFalse(columns & {
            'selected_weight', 'reference_with_weights',
            'weight_ratio', 'weight_resolved'})
