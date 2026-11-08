# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

"""Versioned alignment evidence; no taxonomy or ecological decisions."""

import csv
from fractions import Fraction
import hashlib
import json
import math
import re


SCHEMA_VERSION = 1
QUERY_COLUMNS = ('query_id', 'sequence_length', 'sequence_sha256',
                 'accepted_alignments', 'status')
HIT_COLUMNS = (
    'query_id', 'reference_id', 'strand', 'matches', 'alignment_length',
    'mismatches', 'gap_bases', 'gap_opens', 'gap_extensions', 'query_length',
    'reference_length', 'query_start', 'query_end', 'reference_start',
    'reference_end', 'raw_score', 'cigar', 'identity', 'query_coverage',
    'reference_coverage', 'is_exact_best')
INTEGER_COLUMNS = HIT_COLUMNS[3:16]
SHA256_PATTERN = re.compile(r'^[0-9a-f]{64}$')


def identity_fraction(hit):
    """Exact iddef-2 identity, independent of any printed percentages."""
    return Fraction(int(hit['matches']), int(hit['alignment_length']))


def inventory_fingerprint(rows):
    """Hash ID/length/uppercase-sequence-SHA256 records in supplied order."""
    digest = hashlib.sha256()
    for row in rows:
        digest.update((json.dumps(list(row), ensure_ascii=True,
                                  separators=(',', ':')) + '\n').encode())
    return digest.hexdigest()


def read_tsv(path, columns):
    with open(path, encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames != list(columns):
            raise ValueError('Unexpected columns in %s.' % path)
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('Incorrect number of fields in %s.' % path)
            yield row


def validate_query(row):
    if (not row['query_id'] or int(row['sequence_length']) < 1 or
            int(row['accepted_alignments']) < 0 or
            not SHA256_PATTERN.fullmatch(row['sequence_sha256'])):
        raise ValueError('Invalid query inventory record.')
    expected = 'matched' if int(row['accepted_alignments']) else 'no_hit'
    if row['status'] != expected:
        raise ValueError('Query status disagrees with accepted alignments.')


def validate_hit(row):
    numbers = {key: int(row[key]) for key in INTEGER_COLUMNS}
    if (not row['query_id'] or not row['reference_id'] or
            row['strand'] not in ('+', '-') or
            any(value < 0 for key, value in numbers.items()
                if key != 'raw_score') or
            min(numbers['alignment_length'], numbers['query_length'],
                numbers['reference_length']) < 1):
        raise ValueError('Invalid alignment record.')
    if (numbers['matches'] + numbers['mismatches'] + numbers['gap_bases'] !=
            numbers['alignment_length'] or
            numbers['gap_opens'] + numbers['gap_extensions'] !=
            numbers['gap_bases']):
        raise ValueError('Inconsistent integer alignment statistics.')
    for prefix, length_key in [('query', 'query_length'),
                               ('reference', 'reference_length')]:
        if not (1 <= numbers[prefix + '_start'] <=
                numbers[prefix + '_end'] <= numbers[length_key]):
            raise ValueError('Invalid alignment coordinates.')
    if not re.fullmatch(r'(?:[1-9][0-9]*[MID]|[MID])+', row['cigar']):
        raise ValueError('Invalid VSEARCH CIGAR.')
    operations = [(int(count or 1), operation) for count, operation in
                  re.findall(r'([0-9]*)([MID])', row['cigar'])]
    paired = sum(count for count, operation in operations if operation == 'M')
    query_length = paired + sum(count for count, operation in operations
                                if operation == 'D')
    reference_length = paired + sum(count for count, operation in operations
                                    if operation == 'I')
    internal = operations[:]
    if internal and internal[0][1] != 'M':
        internal.pop(0)
    if internal and internal[-1][1] != 'M':
        internal.pop()
    if (query_length != numbers['query_length'] or
            reference_length != numbers['reference_length'] or
            paired != numbers['matches'] + numbers['mismatches'] or
            sum(count for count, _ in internal) !=
            numbers['alignment_length']):
        raise ValueError('CIGAR disagrees with integer alignment statistics.')
    pairs = numbers['matches'] + numbers['mismatches']
    if pairs > min(numbers['query_length'], numbers['reference_length']):
        raise ValueError('Alignment has more pairs than sequence bases.')
    expected = [float(identity_fraction(row)),
                pairs / numbers['query_length'],
                pairs / numbers['reference_length']]
    for key, value in zip(('identity', 'query_coverage',
                           'reference_coverage'), expected):
        if (not math.isfinite(float(row[key])) or
                abs(float(row[key]) - value) > 1e-12):
            raise ValueError('Derived statistic %s is inconsistent.' % key)
    if row['is_exact_best'] not in ('0', '1'):
        raise ValueError('Invalid exact-best flag.')


def validate_manifest(manifest):
    if manifest.get('schema_version') != SCHEMA_VERSION:
        raise ValueError('Unsupported reference alignment evidence schema.')
    if (manifest.get('identity_definition') != 2 or
            not 0 < manifest.get('identity_floor', 0) <= 1 or
            manifest.get('strand') not in ('plus', 'both') or
            manifest.get('search_mode') not in ('exhaustive', 'heuristic')):
        raise ValueError('Invalid alignment model or search settings.')
    exhaustive = manifest['search_mode'] == 'exhaustive'
    if manifest.get('search_complete') is not exhaustive:
        raise ValueError('Search completeness disagrees with search mode.')
    params = manifest.get('effective_parameters', {})
    required = {'iddef': 2, 'qmask': 'none', 'dbmask': 'dust',
                'minseqlength': 1, 'match': 2, 'mismatch': -4,
                'gapopen': '20I/2E', 'gapext': '2I/1E',
                'wordlength': 8, 'maxhits': 0}
    required.update({'maxaccepts': 0, 'maxrejects': 0, 'minwordmatches': 0}
                    if exhaustive else
                    {'maxaccepts': 1, 'maxrejects': 32,
                     'minwordmatches': 12})
    if any(params.get(key) != value for key, value in required.items()):
        raise ValueError('Effective parameters do not support completeness.')
    if (params.get('id') != manifest['identity_floor'] or
            params.get('strand') != manifest['strand'] or
            not isinstance(params.get('threads'), int) or
            not 0 <= params['threads'] <= 1024 or
            params.get('query_cov') != 0 or params.get('target_cov') != 0):
        raise ValueError('Inconsistent effective search parameters.')
    for key in ('queries', 'references'):
        record = manifest.get(key, {})
        if (record.get('count', 0) < 1 or
                record.get('max_length', 0) < 1 or
                params.get('maxseqlength', 0) < record['max_length'] or
                not SHA256_PATTERN.fullmatch(record.get('sha256', '')) or
                not SHA256_PATTERN.fullmatch(
                    record.get('order_sha256', ''))):
            raise ValueError('Missing sequence inventory fingerprint.')
    if not re.match(r'^vsearch v[0-9]+\.[0-9]+',
                    manifest.get('vsearch_version', '')):
        raise ValueError('Missing VSEARCH version.')
