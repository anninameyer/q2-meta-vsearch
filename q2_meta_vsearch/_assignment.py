# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

"""Unweighted assignment core, independent of QIIME 2 and executables."""

from contextlib import closing
from dataclasses import dataclass
from fractions import Fraction
from itertools import groupby
import json
from pathlib import Path
import sqlite3
import tempfile

from ._alignment_evidence import (
    HIT_COLUMNS, QUERY_COLUMNS, identity_fraction, inventory_fingerprint,
    read_tsv, validate_hit, validate_manifest, validate_query)
from ._taxonomy import consensus_taxonomy, parse_lineage


SUMMARY_COLUMNS = (
    'matched', 'best_identity', 'best_identity_numerator',
    'best_identity_denominator', 'n_accepted_reference_hits',
    'n_best_reference_hits', 'n_best_alignments', 'n_near_best_reference_hits',
    'best_query_coverage', 'best_reference_coverage', 'best_alignment_length',
    'best_mismatches', 'best_gaps', 'best_strand', 'reference_without_weights',
    'n_remaining_reference_hits', 'reference_assignment_status', 'n_best_taxa',
    'taxonomy_without_weights', 'taxonomy_lca', 'taxonomy_ambiguity_label',
    'taxonomic_assignment_status', 'taxonomic_resolution_rank', 'search_mode',
    'search_complete', 'identity_threshold', 'identity_tolerance',
    'reference_selection_policy', 'statistics_reference',
    'n_missing_taxonomy_reference_hits')
NUMERIC_COLUMNS = {
    'best_identity', 'best_identity_numerator', 'best_identity_denominator',
    'n_accepted_reference_hits', 'n_best_reference_hits', 'n_best_alignments',
    'n_near_best_reference_hits', 'best_query_coverage',
    'best_reference_coverage',
    'best_alignment_length', 'best_mismatches', 'best_gaps',
    'n_remaining_reference_hits', 'n_best_taxa', 'identity_threshold',
    'identity_tolerance', 'n_missing_taxonomy_reference_hits'}


@dataclass
class AssignmentResult:
    summary: list
    taxonomy: dict
    mapping: dict


def _parameters(manifest, identity_threshold, identity_tolerance,
                unresolved_policy, max_explicit_species):
    validate_manifest(manifest)
    threshold = Fraction(str(identity_threshold))
    tolerance = Fraction(str(identity_tolerance)) / 100
    if not Fraction(str(manifest['identity_floor'])) <= threshold <= 1:
        raise ValueError('identity_threshold must be at least the stored '
                         'identity floor and at most 1.')
    if not 0 <= tolerance <= 1:
        raise ValueError('identity_tolerance must be 0–100 percentage points.')
    if unresolved_policy not in ('deterministic-reference', 'exclude'):
        raise ValueError('Invalid unresolved_policy.')
    if (not isinstance(max_explicit_species, int) or
            isinstance(max_explicit_species, bool) or
            max_explicit_species < 1):
        raise ValueError('max_explicit_species must be a positive integer.')
    return threshold, tolerance


def resolve_query(fetch_hits, taxonomy, manifest, identity_threshold=0.90,
                  identity_tolerance=0.0,
                  unresolved_policy='deterministic-reference',
                  max_explicit_species=2):
    """Resolve one query from a re-iterable, reference-ID-sorted hit source.

    Only exact-best references enter consensus or reference selection.
    fetch_hits is a callable yielding ordinary dictionaries; no QIIME 2 objects
    or abundance data are required. It is traversed twice without storing all
    alignment rows in memory.
    """
    threshold, tolerance = _parameters(
        manifest, identity_threshold, identity_tolerance, unresolved_policy,
        max_explicit_species)
    maximum = max((identity_fraction(hit) for hit in fetch_hits()),
                  default=Fraction(0))
    row = dict.fromkeys(SUMMARY_COLUMNS)
    row.update(matched='false', n_accepted_reference_hits=0,
               n_best_reference_hits=0, n_best_alignments=0,
               n_near_best_reference_hits=0, n_remaining_reference_hits=0,
               n_best_taxa=0, n_missing_taxonomy_reference_hits=0,
               reference_assignment_status='no_match',
               taxonomy_without_weights='Unassigned',
               taxonomy_lca='Unassigned',
               taxonomy_ambiguity_label='Unassigned',
               taxonomic_assignment_status='no_match',
               taxonomic_resolution_rank='Unassigned',
               search_mode=manifest['search_mode'],
               search_complete=str(manifest['search_complete']).lower(),
               identity_threshold=float(threshold),
               identity_tolerance=float(tolerance * 100),
               reference_selection_policy=unresolved_policy)
    if maximum < threshold:
        return row, 'Unassigned', None
    lineages = set()
    representative = None
    for reference_id, grouped in groupby(
            fetch_hits(), key=lambda hit: hit['reference_id']):
        accepted = False
        exact = []  # At most one hit per strand for one reference.
        near = False
        for hit in grouped:
            identity = identity_fraction(hit)
            accepted |= identity >= threshold
            if identity == maximum:
                exact.append(hit)
            elif identity >= max(threshold, maximum - tolerance):
                near = True
        row['n_accepted_reference_hits'] += int(accepted)
        if not exact:
            row['n_near_best_reference_hits'] += int(near)
            continue
        row['n_best_reference_hits'] += 1
        row['n_best_alignments'] += len(exact)
        if representative is None:
            representative = min(exact, key=lambda hit: hit['strand'] != '+')
        lineage = parse_lineage(taxonomy.get(reference_id))
        if not lineage.values:
            row['n_missing_taxonomy_reference_hits'] += 1
        else:
            lineages.add(lineage)
    lca, label, rank, status, n_taxa = consensus_taxonomy(
        lineages, row['n_missing_taxonomy_reference_hits'],
        max_explicit_species)
    selected = representative['reference_id']
    excluded = (unresolved_policy == 'exclude' and
                status in ('ambiguous', 'missing_taxonomy',
                           'partial_missing_taxonomy'))
    pairs = int(representative['matches']) + int(representative['mismatches'])
    row.update(
        matched='true', best_identity=float(maximum),
        best_identity_numerator=int(representative['matches']),
        best_identity_denominator=int(representative['alignment_length']),
        best_query_coverage=pairs / int(representative['query_length']),
        best_reference_coverage=(
            pairs / int(representative['reference_length'])),
        best_alignment_length=int(representative['alignment_length']),
        best_mismatches=int(representative['mismatches']),
        best_gaps=int(representative['gap_bases']),
        best_strand=representative['strand'], statistics_reference=selected,
        reference_without_weights=None if excluded else selected,
        n_remaining_reference_hits=row['n_best_reference_hits'],
        reference_assignment_status=(
            'excluded_taxonomic_ambiguity' if excluded else
            'unique_reference' if row['n_best_reference_hits'] == 1 else
            'deterministic_reference'),
        n_best_taxa=n_taxa, taxonomy_without_weights=lca, taxonomy_lca=lca,
        taxonomy_ambiguity_label=label, taxonomic_assignment_status=status,
        taxonomic_resolution_rank=rank)
    return row, lca, None if excluded else selected


def assign_reference_alignments(evidence, reference_taxonomy,
                                identity_threshold=0.90,
                                identity_tolerance=0.0,
                                unresolved_policy='deterministic-reference',
                                max_explicit_species=2):
    """Read schema-1 evidence and return summaries, taxonomy, and a mapping.

    Input identifiers are strings; absent taxonomy is reported, never replaced
    by taxonomy from lower-identity hits. Extra taxonomy IDs are allowed:
    many database references need not have accepted alignments.
    """
    root = Path(evidence)
    manifest = json.loads((root / 'search.json').read_text(encoding='utf-8'))
    _parameters(manifest, identity_threshold, identity_tolerance,
                unresolved_policy, max_explicit_species)
    if any(not isinstance(key, str) or not key for key in reference_taxonomy):
        raise ValueError('Reference taxonomy identifiers must be strings.')
    summary, taxonomies, mapping = [], {}, {}
    with tempfile.TemporaryDirectory(prefix='meta-vsearch-assign-') as td:
        with closing(sqlite3.connect(Path(td) / 'hits.sqlite')) as conn:
            conn.execute('CREATE TABLE queries (id TEXT PRIMARY KEY, '
                         'length INT, '
                         'sha256 TEXT, expected INT, observed INT DEFAULT 0)')
            for query in read_tsv(root / 'queries.tsv', QUERY_COLUMNS):
                validate_query(query)
                try:
                    conn.execute('INSERT INTO queries '
                                 '(id, length, sha256, expected) '
                                 'VALUES (?, ?, ?, ?)',
                                 (query['query_id'], query['sequence_length'],
                                  query['sequence_sha256'],
                                  query['accepted_alignments']))
                except sqlite3.IntegrityError as exc:
                    raise ValueError('Duplicate query inventory ID.') from exc
            fingerprint = inventory_fingerprint(conn.execute(
                'SELECT id, length, sha256 FROM queries ORDER BY id'))
            if (fingerprint != manifest['queries']['sha256'] or
                    conn.execute('SELECT count(*) FROM queries').fetchone()[0]
                    != manifest['queries']['count']):
                raise ValueError('Query inventory fingerprint mismatch.')
            conn.execute('CREATE TABLE hits (query_id TEXT, '
                         'reference_id TEXT, '
                         'strand TEXT, data TEXT, '
                         'PRIMARY KEY(query_id, reference_id, strand))')
            for hit in read_tsv(root / 'hits.tsv', HIT_COLUMNS):
                validate_hit(hit)
                query = conn.execute('SELECT length FROM queries WHERE id=?',
                                     (hit['query_id'],)).fetchone()
                if query is None or query[0] != int(hit['query_length']):
                    raise ValueError('Hit does not match query inventory.')
                conn.execute('INSERT INTO hits VALUES (?, ?, ?, ?)',
                             (hit['query_id'], hit['reference_id'],
                              hit['strand'], json.dumps(hit)))
                conn.execute('UPDATE queries SET observed=observed+1 '
                             'WHERE id=?',
                             (hit['query_id'],))
            if conn.execute('SELECT 1 FROM queries WHERE observed != expected '
                            'LIMIT 1').fetchone():
                raise ValueError('Hit counts disagree with query inventory.')
            conn.commit()
            query_ids = conn.execute('SELECT id FROM queries ORDER BY id')
            for query_id, in query_ids:
                def fetch_hits():
                    for data, in conn.execute(
                            'SELECT data FROM hits WHERE query_id=? '
                            'ORDER BY reference_id, strand', (query_id,)):
                        yield json.loads(data)

                row, taxon, reference = resolve_query(
                    fetch_hits, reference_taxonomy, manifest,
                    identity_threshold, identity_tolerance, unresolved_policy,
                    max_explicit_species)
                row = {'id': query_id, **row}
                summary.append(row)
                taxonomies[query_id] = taxon
                if reference is not None:
                    mapping[query_id] = reference
    return AssignmentResult(summary, taxonomies, mapping)
