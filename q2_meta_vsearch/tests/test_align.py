# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

import csv
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from unittest.mock import patch

from qiime2 import Artifact
from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from q2_meta_vsearch._align import align_closed_reference, _search_parameters
from q2_meta_vsearch._alignment_evidence import identity_fraction
from q2_meta_vsearch._format import MetaVsearchReferenceAlignmentsDirFmt


class AlignmentTests(TestPluginBase):
    """Synthetic integration tests against the environment's real VSEARCH."""

    package = 'q2_meta_vsearch.tests'

    def _fasta(self, name):
        return DNAFASTAFormat(
            self.get_data_path('alignment/' + name), mode='r')

    def _run(self, queries='queries.fasta', refs='references.fasta', **kwargs):
        evidence, blast6 = align_closed_reference(
            self._fasta(queries), self._fasta(refs), **kwargs)
        evidence.validate('max')
        blast6.validate()
        root = Path(str(evidence))
        with open(root / 'hits.tsv') as fh:
            hits = list(csv.DictReader(fh, delimiter='\t'))
        with open(root / 'queries.tsv') as fh:
            inventory = list(csv.DictReader(fh, delimiter='\t'))
        manifest = json.loads((root / 'search.json').read_text())
        return evidence, blast6, hits, inventory, manifest

    def test_unique_exact(self):
        _, _, hits, inventory, _ = self._run(
            'unique.fasta', 'reference-single.fasta')
        self.assertEqual(len(hits), 1)
        self.assertEqual(identity_fraction(hits[0]), Fraction(1))
        self.assertEqual(hits[0]['is_exact_best'], '1')
        self.assertEqual(inventory[0]['status'], 'matched')

    def test_all_ties_and_lower_identity_are_retained(self):
        _, _, hits, inventory, _ = self._run()
        best = {hit['reference_id'] for hit in hits
                if hit['is_exact_best'] == '1'}
        self.assertEqual(best, {'ref-a', 'ref-z'})
        self.assertEqual(len(hits), 3)
        lower = next(hit for hit in hits if hit['reference_id'] == 'ref-lower')
        self.assertEqual(identity_fraction(lower), Fraction(99, 100))
        self.assertEqual(lower['is_exact_best'], '0')
        self.assertEqual({row['query_id'] for row in inventory},
                         {'00123', 'unmatched'})

    def test_reference_order_invariance(self):
        _, _, hits, _, first = self._run()
        _, _, reversed_hits, _, second = self._run(
            refs='references-reversed.fasta')
        self.assertEqual(hits, reversed_hits)
        self.assertEqual(first['references']['sha256'],
                         second['references']['sha256'])
        self.assertNotEqual(first['references']['order_sha256'],
                            second['references']['order_sha256'])

    def test_first_acceptable_can_be_suboptimal(self):
        _, _, exhaustive, _, _ = self._run(
            'suboptimal-query.fasta', 'suboptimal-references.fasta')
        _, _, heuristic, _, manifest = self._run(
            'suboptimal-query.fasta', 'suboptimal-references.fasta',
            search_mode='heuristic')
        self.assertEqual(len(exhaustive), 2)
        self.assertEqual(len(heuristic), 1)
        self.assertEqual(heuristic[0]['reference_id'], 'worse')
        self.assertLess(identity_fraction(heuristic[0]),
                        max(map(identity_fraction, exhaustive)))
        self.assertFalse(manifest['search_complete'])

    def test_exhaustive_agrees_with_single_reference_oracle(self):
        _, _, combined, _, _ = self._run(
            'suboptimal-query.fasta', 'suboptimal-references.fasta')
        individually = []
        for name in ['suboptimal-better.fasta', 'suboptimal-worse.fasta']:
            _, _, hits, _, _ = self._run('suboptimal-query.fasta', name)
            individually.extend(hits)
        oracle = {hit['reference_id']: identity_fraction(hit)
                  for hit in individually}
        self.assertEqual({hit['reference_id']: identity_fraction(hit)
                          for hit in combined}, oracle)
        self.assertEqual(max(oracle.values()), Fraction(99, 100))
        self.assertEqual({hit['reference_id'] for hit in combined
                          if hit['is_exact_best'] == '1'}, {'better'})

    def test_exhaustive_bypasses_word_filter(self):
        _, _, exhaustive, _, _ = self._run(
            'unique.fasta', 'reference-ambiguous.fasta')
        _, _, heuristic, _, _ = self._run(
            'unique.fasta', 'reference-ambiguous.fasta',
            search_mode='heuristic')
        self.assertEqual(len(exhaustive), 1)
        self.assertEqual(identity_fraction(exhaustive[0]), Fraction(1))
        self.assertEqual(heuristic, [])

    def test_exact_precision_hidden_by_blast6(self):
        _, blast6, hits, _, _ = self._run(
            'precision-query.fasta', 'precision-references.fasta')
        with open(str(blast6)) as fh:
            rounded = [row[2] for row in csv.reader(fh, delimiter='\t')]
        self.assertEqual(rounded, ['99.9', '99.9'])
        self.assertNotEqual(identity_fraction(hits[0]),
                            identity_fraction(hits[1]))
        self.assertEqual([hit['reference_id'] for hit in hits
                          if hit['is_exact_best'] == '1'], ['better'])

    def test_strands_gaps_coverage_and_ambiguous_symbols(self):
        _, _, hits, inventory, _ = self._run(
            'geometry.fasta', 'reference-single.fasta', strand='both')
        by_query = {hit['query_id']: hit for hit in hits}
        self.assertEqual(len(by_query), 5)
        self.assertEqual(by_query['forward']['strand'], '+')
        self.assertEqual(by_query['reverse']['strand'], '-')
        self.assertEqual(identity_fraction(by_query['ambiguous']), Fraction(1))
        short = by_query['short']
        self.assertEqual(short['alignment_length'], '160')
        self.assertEqual(float(short['query_coverage']), 1)
        self.assertEqual(float(short['reference_coverage']), 0.8)
        self.assertEqual(short['gap_bases'], '0')  # excludes terminal gaps
        self.assertIn('I', short['cigar'])
        gap = by_query['internal-gap']
        self.assertEqual(int(gap['gap_bases']), 14)
        self.assertEqual(int(gap['gap_opens']), 1)
        self.assertEqual(int(gap['gap_extensions']), 13)
        self.assertEqual(gap['alignment_length'], '214')
        self.assertEqual(len(inventory), 5)

    def test_reverse_does_not_match_on_plus(self):
        _, _, hits, inventory, _ = self._run(
            'geometry.fasta', 'reference-single.fasta')
        self.assertNotIn('reverse', {hit['query_id'] for hit in hits})
        reverse = next(row for row in inventory
                       if row['query_id'] == 'reverse')
        self.assertEqual(reverse['status'], 'no_hit')

    def test_no_hits_have_inventory_and_no_fake_alignments(self):
        _, blast6, hits, inventory, _ = self._run(
            'unmatched.fasta', 'reference-single.fasta')
        self.assertEqual(hits, [])
        self.assertEqual(inventory[0]['accepted_alignments'], '0')
        self.assertEqual(inventory[0]['status'], 'no_hit')
        with open(str(blast6)) as fh:
            self.assertEqual(fh.read().strip().split('\t'),
                             ['unmatched', '*', '0.0', '0', '0', '0',
                              '0', '0', '0', '0', '-1', '0'])

    def test_problematic_identifiers_preserved_exactly(self):
        _, blast6, hits, inventory, _ = self._run(
            'identifiers.fasta', 'identifiers.fasta')
        expected = {'00123', 'NA', '1e3', 'asv;size=7;odd', 'asv|x:y'}
        self.assertEqual({row['query_id'] for row in inventory}, expected)
        self.assertEqual({hit['query_id'] for hit in hits}, expected)
        self.assertEqual({hit['reference_id'] for hit in hits}, expected)
        with open(str(blast6)) as fh:
            rows = list(csv.reader(fh, delimiter='\t'))
        self.assertEqual({row[0] for row in rows}, expected)
        self.assertEqual({row[1] for row in rows}, expected)

    def test_identity_floor(self):
        _, _, hits, _, manifest = self._run(identity_floor=1.0)
        self.assertEqual(len(hits), 2)
        self.assertTrue(all(identity_fraction(hit) == 1 for hit in hits))
        self.assertEqual(manifest['identity_floor'], 1.0)
        _, _, boundary_hits, _, _ = self._run(identity_floor=0.99)
        self.assertEqual(len(boundary_hits), 3)

    def test_manifest_and_fingerprints(self):
        _, _, _, inventory, manifest = self._run(threads=2)
        self.assertEqual(manifest['schema_version'], 1)
        self.assertEqual(manifest['identity_definition'], 2)
        self.assertEqual(manifest['identity_floor'], 0.9)
        self.assertTrue(manifest['search_complete'])
        self.assertEqual(manifest['search_mode'], 'exhaustive')
        params = manifest['effective_parameters']
        for key in ('maxaccepts', 'maxrejects', 'minwordmatches', 'maxhits'):
            self.assertEqual(params[key], 0)
        self.assertEqual(params['threads'], 2)
        self.assertEqual(params['qmask'], 'none')
        self.assertEqual(params['dbmask'], 'dust')
        self.assertEqual(manifest['references']['count'], 3)
        self.assertRegex(manifest['vsearch_version'], r'^vsearch v\d+\.')
        self.assertGreater(manifest['vsearch_elapsed_seconds'], 0)
        source = Path(self.get_data_path('alignment/queries.fasta'))
        sequence = source.read_text().splitlines()[1]
        self.assertEqual(inventory[0]['sequence_sha256'],
                         hashlib.sha256(sequence.encode()).hexdigest())

    def test_registered_action_and_artifact_roundtrip(self):
        action = self.plugin.methods['align_closed_reference']
        sequences = Artifact.import_data('FeatureData[Sequence]',
                                         self._fasta('queries.fasta'))
        refs = Artifact.import_data('FeatureData[Sequence]',
                                    self._fasta('references.fasta'))
        outputs = action(sequences=sequences, reference_sequences=refs)
        self.assertEqual(str(outputs.alignments.type),
                         'FeatureData[ReferenceAlignments]')
        self.assertEqual(str(outputs.blast6.type), 'FeatureData[BLAST6]')
        outputs.alignments.validate(level='max')
        path = outputs.alignments.save(str(Path(self.temp_dir.name) / 'hits'))
        loaded = Artifact.load(path)
        loaded.validate(level='max')
        view = loaded.view(MetaVsearchReferenceAlignmentsDirFmt)
        self.assertTrue((Path(str(view)) / 'search.json').is_file())

    def test_all_unmatched_registered_action(self):
        method = self.plugin.methods['align_closed_reference']
        outputs = method(
            sequences=Artifact.import_data(
                'FeatureData[Sequence]', self._fasta('unmatched.fasta')),
            reference_sequences=Artifact.import_data(
                'FeatureData[Sequence]',
                self._fasta('reference-single.fasta')))
        outputs.alignments.validate(level='max')
        outputs.blast6.validate(level='max')

    def test_corrupt_evidence_rejected(self):
        evidence, _, _, _, _ = self._run()
        destination = Path(self.temp_dir.name) / 'corrupt'
        shutil.copytree(str(evidence), destination)
        hits = destination / 'hits.tsv'
        hits.write_text(hits.read_text().replace('\t1\n', '\t0\n'))
        with self.assertRaisesRegex(ValidationError, 'exact-best'):
            MetaVsearchReferenceAlignmentsDirFmt(
                destination, mode='r').validate('max')

    def test_manifest_contradictions_rejected(self):
        evidence, _, _, _, manifest = self._run()
        changes = [
            ('schema_version', 999), ('search_complete', False),
            ('identity_definition', 1), ('identity_floor', 1.1),
            ('queries', {**manifest['queries'], 'sha256': '0' * 64}),
            ('effective_parameters',
             {**manifest['effective_parameters'], 'maxaccepts': 1})]
        for index, (key, value) in enumerate(changes):
            with self.subTest(field=key):
                destination = Path(self.temp_dir.name) / str(index)
                shutil.copytree(str(evidence), destination)
                altered = copy.deepcopy(manifest)
                altered[key] = value
                (destination / 'search.json').write_text(json.dumps(altered))
                with self.assertRaises(ValidationError):
                    MetaVsearchReferenceAlignmentsDirFmt(
                        destination, mode='r').validate('max')

    def test_hit_statistics_and_membership_rejected(self):
        evidence, _, hits, _, _ = self._run()
        changes = [('matches', '0'), ('query_coverage', 'nan'),
                   ('query_id', 'unknown'), ('cigar', '200I')]
        for index, (key, value) in enumerate(changes):
            with self.subTest(field=key):
                destination = Path(self.temp_dir.name) / str(index)
                shutil.copytree(str(evidence), destination)
                altered = copy.deepcopy(hits)
                altered[0][key] = value
                with open(destination / 'hits.tsv', 'w') as fh:
                    writer = csv.DictWriter(fh, fieldnames=altered[0],
                                            delimiter='\t')
                    writer.writeheader()
                    writer.writerows(altered)
                with self.assertRaises(ValidationError):
                    MetaVsearchReferenceAlignmentsDirFmt(
                        destination, mode='r').validate('max')

    def test_query_counts_and_duplicate_hits_rejected(self):
        evidence, _, _, _, _ = self._run()
        for index, name in enumerate(['queries.tsv', 'hits.tsv']):
            destination = Path(self.temp_dir.name) / str(index)
            shutil.copytree(str(evidence), destination)
            path = destination / name
            lines = path.read_text().splitlines(keepends=True)
            path.write_text(''.join(lines + [lines[1]]))
            with self.assertRaises(ValidationError):
                MetaVsearchReferenceAlignmentsDirFmt(
                    destination, mode='r').validate('max')

    def test_duplicate_and_empty_inputs_rejected(self):
        for name, message in [('duplicate-ids.fasta', 'Duplicate'),
                              ('empty.fasta', 'must not be empty')]:
            with self.assertRaisesRegex(ValueError, message):
                align_closed_reference(self._fasta(name),
                                       self._fasta('references.fasta'))

    def test_external_failure_is_clear(self):
        with patch('q2_meta_vsearch._align.run_command') as run:
            run.side_effect = subprocess.CalledProcessError(7, ['vsearch'])
            with self.assertRaisesRegex(RuntimeError, 'exit code 7'):
                align_closed_reference(self._fasta('unique.fasta'),
                                       self._fasta('references.fasta'))

    def test_invalid_parameters_rejected(self):
        base = dict(identity_floor=0.9, strand='plus',
                    search_mode='exhaustive', threads=1, maximum=200)
        for key, value in [('identity_floor', 0), ('identity_floor', 1.1),
                           ('strand', 'minus'), ('search_mode', 'unknown'),
                           ('threads', -1)]:
            with self.subTest(parameter=key, value=value):
                with self.assertRaises(ValueError):
                    _search_parameters(**{**base, key: value})
