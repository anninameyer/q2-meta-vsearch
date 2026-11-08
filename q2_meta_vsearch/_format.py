# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv
from contextlib import closing
import json
from pathlib import Path
import sqlite3
import tempfile

import numpy as np
import qiime2.plugin.model as model
from qiime2.plugin import ValidationError

from ._alignment_evidence import (
    HIT_COLUMNS, QUERY_COLUMNS, identity_fraction, inventory_fingerprint,
    read_tsv, validate_hit, validate_manifest, validate_query)


class MetaVsearchUchimeStatsFmt(model.TextFileFormat):

    def _check_n_records(self, n):
        with open(str(self)) as fh:
            csv_reader = csv.reader(fh, delimiter='\t')
            for i, row in enumerate(csv_reader):
                if i == n:
                    break
                else:
                    if len(row) != 18:
                        raise ValidationError(
                            'Incorrect number of fields detected on line %d.'
                            ' Should be exactly 18.' % (i + 1))

    def _validate_(self, level):
        record_count_map = {'min': 5, 'max': np.inf}
        self._check_n_records(record_count_map[level])


MetaVsearchUchimeStatsDirFmt = model.SingleFileDirectoryFormat(
    'MetaVsearchUchimeStatsDirFmt', 'stats.tsv', MetaVsearchUchimeStatsFmt)


class MetaVsearchReferenceHitsFormat(model.TextFileFormat):
    def _validate_(self, level):
        try:
            for index, row in enumerate(read_tsv(str(self), HIT_COLUMNS)):
                validate_hit(row)
                if level == 'min' and index >= 4:
                    break
        except (ValueError, KeyError, TypeError, ArithmeticError) as exc:
            raise ValidationError(str(exc)) from exc


class MetaVsearchQueryInventoryFormat(model.TextFileFormat):
    def _validate_(self, level):
        try:
            count = 0
            for count, row in enumerate(read_tsv(str(self), QUERY_COLUMNS), 1):
                validate_query(row)
                if level == 'min' and count >= 5:
                    break
            if not count:
                raise ValueError('Query inventory must not be empty.')
        except (ValueError, KeyError, TypeError) as exc:
            raise ValidationError(str(exc)) from exc


class MetaVsearchSearchManifestFormat(model.TextFileFormat):
    def _validate_(self, level):
        try:
            with self.open() as fh:
                validate_manifest(json.load(fh))
        except (ValueError, KeyError, TypeError, AttributeError) as exc:
            raise ValidationError(str(exc)) from exc


class MetaVsearchReferenceAlignmentsDirFmt(model.DirectoryFormat):
    hits = model.File('hits.tsv', format=MetaVsearchReferenceHitsFormat)
    queries = model.File('queries.tsv', format=MetaVsearchQueryInventoryFormat)
    search = model.File('search.json', format=MetaVsearchSearchManifestFormat)

    def _validate_(self, level):
        if level == 'min':
            return
        # Cross-file validation uses disk space rather than memory proportional
        # to the (potentially very large) number of accepted alignments.
        root = Path(str(self))
        with tempfile.TemporaryDirectory(prefix='meta-vsearch-check-') as td:
            with closing(sqlite3.connect(Path(td) / 'check.sqlite')) as conn:
                try:
                    self._check_consistency(root, conn)
                except (ValueError, KeyError, TypeError, ArithmeticError,
                        sqlite3.IntegrityError) as exc:
                    raise ValidationError(str(exc)) from exc

    @staticmethod
    def _check_consistency(root, conn):
        manifest = json.loads((root / 'search.json').read_text())
        conn.execute('CREATE TABLE queries (id TEXT PRIMARY KEY, length INT, '
                     'sha256 TEXT, expected INT, observed INT, best_m INT, '
                     'best_l INT)')
        for row in read_tsv(root / 'queries.tsv', QUERY_COLUMNS):
            conn.execute('INSERT INTO queries VALUES (?, ?, ?, ?, 0, 0, 1)',
                         (row['query_id'], int(row['sequence_length']),
                          row['sequence_sha256'],
                          int(row['accepted_alignments'])))
        count = conn.execute('SELECT count(*) FROM queries').fetchone()[0]
        fingerprint = inventory_fingerprint(conn.execute(
            'SELECT id, length, sha256 FROM queries ORDER BY id'))
        if (count != manifest['queries']['count'] or
                fingerprint != manifest['queries']['sha256']):
            raise ValueError('Query inventory fingerprint does not match.')
        conn.execute('CREATE TABLE pairs (query_id TEXT, reference_id TEXT, '
                     'strand TEXT, '
                     'PRIMARY KEY(query_id, reference_id, strand))')
        for row in read_tsv(root / 'hits.tsv', HIT_COLUMNS):
            query = conn.execute('SELECT length, best_m, best_l FROM queries '
                                 'WHERE id=?', (row['query_id'],)).fetchone()
            if query is None or query[0] != int(row['query_length']):
                raise ValueError('Hit does not match query inventory.')
            if (row['strand'] == '-' and manifest['strand'] != 'both' or
                    float(identity_fraction(row)) + 1e-15 <
                    manifest['identity_floor']):
                raise ValueError('Hit contradicts recorded search settings.')
            conn.execute('INSERT INTO pairs VALUES (?, ?, ?)',
                         (row['query_id'], row['reference_id'], row['strand']))
            matches, length = int(row['matches']), int(row['alignment_length'])
            conn.execute('UPDATE queries SET observed=observed+1 WHERE id=?',
                         (row['query_id'],))
            if matches * query[2] > query[1] * length:
                conn.execute('UPDATE queries SET best_m=?, best_l=? '
                             'WHERE id=?',
                             (matches, length, row['query_id']))
        if conn.execute('SELECT 1 FROM queries WHERE observed != expected '
                        'LIMIT 1').fetchone():
            raise ValueError('Hit counts do not match query inventory.')
        for row in read_tsv(root / 'hits.tsv', HIT_COLUMNS):
            best = conn.execute('SELECT best_m, best_l FROM queries '
                                'WHERE id=?',
                                (row['query_id'],)).fetchone()
            exact = int(row['matches']) * best[1] == \
                best[0] * int(row['alignment_length'])
            if int(row['is_exact_best']) != int(exact):
                raise ValueError('Incorrect exact-best alignment flag.')
