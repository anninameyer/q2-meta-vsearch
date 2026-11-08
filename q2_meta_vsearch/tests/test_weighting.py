# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

import csv
from hashlib import sha256
from pathlib import Path
import platform
from unittest.mock import patch

import biom
import numpy as np
import pandas as pd
from qiime2 import Artifact, Metadata
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import TSVTaxonomyFormat

from q2_meta_vsearch._assign import read_class_weights, read_reference_taxonomy
from q2_meta_vsearch._assignment import (
    assign_reference_alignments, SUMMARY_COLUMNS, UNWEIGHTED_COLUMNS)
from q2_meta_vsearch._format import MetaVsearchReferenceAlignmentsDirFmt
from q2_meta_vsearch._weighting import (
    prepare_weights, WEIGHT_NUMERIC_COLUMNS)


class WeightingTests(TestPluginBase):
    package = 'q2_meta_vsearch.tests'

    def setUp(self):
        super().setUp()
        platform.platform()
        self.external = patch('subprocess.run', side_effect=AssertionError(
            'Weighting must not execute external commands')).start()
        self.addCleanup(patch.stopall)

    def _path(self, name):
        return Path(self.get_data_path('assignment/' + name))

    def _weights(self, scenario='default', fixture='weights.tsv'):
        with self._path(fixture).open() as fh:
            return {r['Taxon']: float(r[scenario])
                    for r in csv.DictReader(fh, delimiter='\t') if r[scenario]}

    def _assign(self, scenario='default', taxonomy='two-species.tsv',
                evidence='evidence', weights=None, **kwargs):
        return assign_reference_alignments(
            self._path(evidence),
            read_reference_taxonomy(self._path(taxonomy)),
            class_weights=(self._weights(scenario) if weights is None
                           else weights),
            **kwargs)

    def _row(self, result, query='00123'):
        return next(r for r in result.summary if r['id'] == query)

    def _table(self, weights=None, sample='Weight'):
        weights = self._weights() if weights is None else weights
        return biom.Table(np.array(list(weights.values()))[:, None],
                          list(weights), [sample])

    def _action(self, table=None, evidence='evidence', **kwargs):
        alignments = Artifact.import_data(
            'FeatureData[ReferenceAlignments]', str(self._path(evidence)),
            view_type=MetaVsearchReferenceAlignmentsDirFmt)
        taxa = Artifact.import_data('FeatureData[Taxonomy]',
                                    str(self._path('two-species.tsv')),
                                    view_type=TSVTaxonomyFormat)
        weights = Artifact.import_data('FeatureTable[RelativeFrequency]',
                                       self._table() if table is None
                                       else table)
        return self.plugin.methods['assign_closed_reference'](
            alignments=alignments, reference_taxonomy=taxa,
            class_weights=weights, **kwargs)

    def test_unique_highest_taxon_and_metadata(self):
        result = self._assign()
        row = self._row(result)
        self.assertEqual(row['reference_without_weights'], 'ref-a')
        self.assertEqual(row['reference_with_weights'], 'ref-c')
        self.assertEqual(result.mapping['00123'], 'ref-c')
        self.assertEqual(row['n_best_reference_hits'], 3)
        self.assertEqual(row['n_best_taxa'], 2)
        self.assertEqual(row['n_weighted_best_reference_hits'], 1)
        self.assertEqual(row['n_weighted_best_taxa'], 1)
        self.assertEqual(row['selected_weight'], .2)
        self.assertEqual(row['second_best_weight'], .05)
        self.assertEqual(row['weight_difference'], .15)
        self.assertEqual(row['weight_ratio'], 4)
        for key in ['weight_resolved', 'assignment_changed',
                    'taxonomy_changed']:
            self.assertEqual(row[key], 'true')
        self.assertNotIn('s__', row['taxonomy_without_weights'])
        self.assertIn('s__Lactiplantibacillus paraplantarum',
                      row['taxonomy_with_weights'])
        self.assertEqual(result.taxonomy['00123'],
                         row['taxonomy_with_weights'])
        self.external.assert_not_called()

    def test_lower_identity_and_near_best_never_eligible(self):
        row = self._row(self._assign(identity_tolerance=1))
        self.assertEqual(row['n_near_best_reference_hits'], 1)
        self.assertEqual(row['reference_with_weights'], 'ref-c')
        self.assertNotIn('Archaea', row['taxonomy_with_weights'])
        for threshold in range(90, 101):
            with self.subTest(threshold=threshold):
                r = self._row(self._assign(identity_threshold=threshold / 100))
                self.assertEqual(r['reference_with_weights'], 'ref-c')

    def test_equal_and_all_zero_weights_retain_all(self):
        for scenario in ['equal', 'zeros']:
            with self.subTest(scenario=scenario):
                row = self._row(self._assign(scenario))
                self.assertEqual(row['reference_with_weights'], 'ref-a')
                self.assertEqual(row['n_weighted_best_reference_hits'], 3)
                self.assertEqual(row['n_weighted_best_taxa'], 2)
                self.assertEqual(row['weight_assignment_status'],
                                 'equal_weights')
                self.assertEqual(row['weight_resolved'], 'false')
                self.assertIsNone(row['selected_weight'])
                self.assertEqual(row['weight_difference'], 0)
                self.assertEqual(row['weight_ratio'], 1)
                self.assertEqual(row['taxonomy_with_weights'],
                                 row['taxonomy_without_weights'])

    def test_zero_runner_and_extremely_small_weights(self):
        for scenario in ['positive-zero', 'tiny', 'overflow']:
            with self.subTest(scenario=scenario):
                row = self._row(self._assign(scenario))
                self.assertEqual(row['reference_with_weights'], 'ref-c')
                self.assertEqual(row['weight_assignment_status'], 'resolved')
                self.assertEqual(row['weight_resolved'], 'true')
                if scenario != 'tiny':
                    self.assertIsNone(row['weight_ratio'])
                else:
                    self.assertEqual(row['weight_ratio'], 2)
        row = self._row(self._assign('positive-zero', min_weight_ratio=1e300))
        self.assertEqual(row['reference_with_weights'], 'ref-c')

    def test_both_margins_and_exact_decimal_boundary(self):
        for difference, ratio, expected in [
                (.6, 4, 'resolved'), (.600001, 4, 'insufficient_margin'),
                (.6, 4.00001, 'insufficient_margin'),
                (0, 5, 'insufficient_margin'), (.7, 1, 'insufficient_margin')]:
            with self.subTest(difference=difference, ratio=ratio):
                row = self._row(self._assign(
                    'margin', min_weight_difference=difference,
                    min_weight_ratio=ratio))
                self.assertEqual(row['weight_assignment_status'], expected)
                if expected == 'insufficient_margin':
                    self.assertEqual(row['n_weighted_best_taxa'], 2)
                    self.assertEqual(row['reference_with_weights'], 'ref-a')
        # Original .2 - .05 fails .2; per-ASV renormalization would pass.
        row = self._row(self._assign(min_weight_difference=.2))
        self.assertEqual(row['weight_assignment_status'],
                         'insufficient_margin')

    def test_missing_weights_default_and_error(self):
        row = self._row(self._assign('missing'))
        self.assertEqual(row['missing_weight_count'], 1)
        self.assertEqual(row['weight_assignment_status'], 'missing_weights')
        self.assertEqual(row['n_weighted_best_reference_hits'], 3)
        self.assertEqual(row['reference_with_weights'], 'ref-a')
        self.assertIsNone(row['selected_weight'])
        self.assertIsNone(row['second_best_weight'])
        with self.assertRaisesRegex(
                ValueError, 'Missing ecological weights.*lineages'):
            self._assign('missing', missing_weight_policy='error')

    def test_missing_taxonomy_cannot_be_weighted(self):
        row = self._row(self._assign(taxonomy='missing-taxonomy.tsv'))
        self.assertEqual(row['weight_assignment_status'], 'missing_taxonomy')
        self.assertGreater(row['missing_weight_count'], 0)
        self.assertGreater(row['n_missing_taxonomy_reference_hits'], 0)
        self.assertEqual(row['taxonomy_with_weights'], 'Unassigned')
        self.assertEqual(row['n_weighted_best_reference_hits'], 3)
        with self.assertRaisesRegex(
                ValueError, 'missing/unusable reference taxonomy'):
            self._assign(taxonomy='missing-taxonomy.tsv',
                         missing_weight_policy='error')

    def test_incompatible_databases_do_not_match_terminal_species(self):
        weights = self._weights('Weight', 'weights-incompatible.tsv')
        row = self._row(self._assign(weights=weights))
        self.assertEqual(row['missing_weight_count'], 2)
        self.assertEqual(row['taxonomy_with_weights'], row['taxonomy_lca'])
        self.assertEqual(row['weight_assignment_status'], 'missing_weights')

    def test_whitespace_normalization_and_collision_detection(self):
        weights = self._weights()
        spaced = {'  ' + key.replace('; ', ' ;   ') + ' ': value
                  for key, value in weights.items()}
        row = self._row(self._assign(weights=spaced))
        self.assertEqual(row['reference_with_weights'], 'ref-c')
        key = next(iter(weights))
        weights[' ' + key] = weights[key]
        with self.assertRaisesRegex(ValueError, 'normalization collision'):
            self._assign(weights=weights)

    def test_reference_multiplicity_and_order_independence(self):
        baseline = self._row(self._assign())
        more = baseline
        baseline = self._row(self._assign(evidence='evidence-less-duplicates'))
        self.assertGreater(more['n_best_reference_hits'],
                           baseline['n_best_reference_hits'])
        for key in ['reference_with_weights', 'taxonomy_with_weights',
                    'selected_weight', 'weight_ratio']:
            self.assertEqual(more[key], baseline[key])
        permuted = self._assign(evidence='evidence-reversed')
        self.assertEqual(permuted.summary, self._assign().summary)

    def test_same_taxon_is_not_new_taxonomic_resolution(self):
        row = self._row(self._assign(taxonomy='same-species.tsv'))
        self.assertEqual(row['n_best_reference_hits'], 3)
        self.assertEqual(row['n_weighted_best_taxa'], 1)
        self.assertEqual(row['n_weighted_best_reference_hits'], 3)
        self.assertEqual(row['weight_assignment_status'], 'already_resolved')
        self.assertEqual(row['weight_resolved'], 'false')
        self.assertEqual(row['taxonomy_changed'], 'false')
        self.assertEqual(row['reference_with_weights'], 'ref-a')
        self.assertEqual(row['selected_weight'], .05)
        self.assertIsNone(row['second_best_weight'])

    def test_unique_reference_and_selection_within_winning_taxon(self):
        row = self._row(self._assign(evidence='evidence-unique'), 'unique')
        self.assertEqual(row['reference_with_weights'], 'ref-a')
        self.assertEqual(row['weight_assignment_status'], 'already_resolved')
        weights = self._weights()
        taxa = read_reference_taxonomy(self._path('two-species.tsv'))
        weights[taxa['ref-a']], weights[taxa['ref-c']] = .2, .05
        row = self._row(self._assign(weights=weights))
        self.assertEqual(row['reference_with_weights'], 'ref-a')
        self.assertEqual(row['n_weighted_best_reference_hits'], 2)
        self.assertEqual(row['n_weighted_best_taxa'], 1)
        self.assertEqual(row['weight_resolved'], 'true')

    def test_missing_lower_identity_weights_do_not_raise(self):
        weights = self._weights()
        taxa = read_reference_taxonomy(self._path('two-species.tsv'))
        del weights[taxa['ref-lower']]
        row = self._row(self._assign(weights=weights,
                                     missing_weight_policy='error',
                                     identity_tolerance=1))
        self.assertEqual(row['missing_weight_count'], 0)
        self.assertEqual(row['reference_with_weights'], 'ref-c')

    def test_incomplete_lineages_do_not_invent_species(self):
        taxa = read_reference_taxonomy(self._path('missing-species.tsv'))
        row = self._row(self._assign(taxonomy='missing-species.tsv',
                                     weights={taxon: .25
                                              for taxon in taxa.values()}))
        self.assertNotIn('s__', row['taxonomy_with_weights'])
        self.assertEqual(row['weight_resolved'], 'false')

    def test_tied_top_retains_all_candidates_conservatively(self):
        taxa = read_reference_taxonomy(self._path('three-species.tsv'))
        weights = {taxa['ref-a']: .4, taxa['ref-b']: .2,
                   taxa['ref-c']: .4, taxa['ref-lower']: 0}
        row = self._row(self._assign(taxonomy='three-species.tsv',
                                     weights=weights))
        self.assertEqual(row['weight_assignment_status'], 'equal_weights')
        self.assertEqual(row['n_weighted_best_taxa'], 3)
        self.assertEqual(row['n_weighted_best_reference_hits'], 3)
        self.assertEqual(row['taxonomy_with_weights'], row['taxonomy_lca'])

    def test_exclusion_uses_final_ambiguity(self):
        resolved = self._row(self._assign(unresolved_policy='exclude'))
        self.assertIsNone(resolved['reference_without_weights'])
        self.assertEqual(resolved['reference_with_weights'], 'ref-c')
        tied = self._row(self._assign('equal', unresolved_policy='exclude'))
        self.assertIsNone(tied['reference_with_weights'])
        self.assertEqual(tied['n_weighted_best_reference_hits'], 3)
        missing = self._row(self._assign(
            'missing', unresolved_policy='exclude'))
        self.assertIsNone(missing['reference_with_weights'])

    def test_unmatched_and_original_ids(self):
        result = self._assign()
        self.assertEqual({r['id'] for r in result.summary},
                         {'00123', 'NA', '1e3', 'asv;size=7;odd',
                          'boundary', 'unmatched'})
        self.assertEqual(len(result.summary), 6)
        row = self._row(result, 'unmatched')
        self.assertEqual(row['weight_assignment_status'], 'no_match')
        self.assertEqual(row['n_weighted_best_taxa'], 0)
        self.assertIsNone(row['reference_with_weights'])
        self.assertIsNone(row['selected_weight'])
        self.assertEqual(row['taxonomy_with_weights'], 'Unassigned')

    def test_unweighted_columns_and_alignment_files_unchanged(self):
        root = self._path('evidence')
        before = {p.name: sha256(p.read_bytes()).hexdigest()
                  for p in root.iterdir()}
        taxa = read_reference_taxonomy(self._path('two-species.tsv'))
        baseline = assign_reference_alignments(root, taxa)
        for scenario in ['default', 'equal', 'missing']:
            weighted = self._assign(scenario)
            for old, new in zip(baseline.summary, weighted.summary):
                self.assertEqual({key: old[key] for key in UNWEIGHTED_COLUMNS},
                                 {key: new[key] for key in UNWEIGHTED_COLUMNS})
        self.assertEqual(before, {p.name: sha256(p.read_bytes()).hexdigest()
                                  for p in root.iterdir()})
        self.assertEqual(assign_reference_alignments(root, taxa), baseline)
        self.external.assert_not_called()

    def test_invalid_weights_and_parameters(self):
        key = next(iter(self._weights()))
        for value in [-1, float('nan'), float('inf'), True]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, 'finite.*nonnegative'):
                    prepare_weights({key: value})
        with self.assertRaisesRegex(ValueError, 'at least one taxon'):
            self._assign(weights={})
        for kwargs in [dict(min_weight_ratio=.9),
                       dict(min_weight_ratio=float('inf')),
                       dict(min_weight_difference=-1),
                       dict(min_weight_difference=float('nan')),
                       dict(missing_weight_policy='zero')]:
            with self.subTest(kwargs=kwargs):
                with self.assertRaises(ValueError):
                    self._assign(**kwargs)

    def test_single_vector_with_arbitrary_sample_identifier(self):
        self.assertEqual(read_class_weights(self._table(sample='Habitat')),
                         self._weights())
        two = self._table().merge(self._table(sample='Other'))
        with self.assertRaisesRegex(ValueError, 'exactly one weight vector'):
            read_class_weights(two)
        with self.assertRaisesRegex(ValueError, 'exactly one weight vector'):
            self._action(table=two)

    def test_weighted_action_artifact_roundtrip(self):
        outputs = self._action()
        outputs.assignment_metadata.validate(level='max')
        outputs.taxonomy.validate(level='max')
        md = outputs.assignment_metadata.view(Metadata)
        frame = md.to_dataframe()
        self.assertEqual(list(frame.columns),
                         list(SUMMARY_COLUMNS))
        for column in WEIGHT_NUMERIC_COLUMNS:
            self.assertEqual(md.columns[column].type, 'numeric')
        self.assertEqual(frame.loc['00123', 'reference_with_weights'], 'ref-c')
        self.assertEqual(frame.loc['00123', 'reference_without_weights'],
                         'ref-a')
        taxa = outputs.taxonomy.view(pd.DataFrame)
        self.assertEqual(taxa.loc['00123', 'Taxon'],
                         frame.loc['00123', 'taxonomy_with_weights'])
        self.assertEqual(set(frame.index), set(taxa.index))
        path = str(Path(self.temp_dir.name) / 'weighted.tsv')
        md.save(path)
        pd.testing.assert_frame_equal(
            Metadata.load(path).to_dataframe(), frame)
        saved = outputs.assignment_metadata.save(
            str(Path(self.temp_dir.name) / 'weighted'))
        pd.testing.assert_frame_equal(
            Artifact.load(saved).view(pd.DataFrame), frame)
        self.external.assert_not_called()

    def test_metadata_handles_infinite_ratio_and_all_unmatched(self):
        for scenario in ['positive-zero', 'overflow']:
            outputs = self._action(table=self._table(self._weights(scenario)))
            outputs.assignment_metadata.validate(level='max')
            self.assertTrue(pd.isna(outputs.assignment_metadata.view(
                pd.DataFrame).loc['00123', 'weight_ratio']))
        outputs = self._action(evidence='evidence-unmatched')
        outputs.assignment_metadata.validate(level='max')
        frame = outputs.assignment_metadata.view(pd.DataFrame)
        self.assertTrue(frame['selected_weight'].isna().all())
