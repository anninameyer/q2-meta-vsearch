# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

"""Table-independent VSEARCH alignment with reusable sequence evidence."""

import csv
from contextlib import closing
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import subprocess
import tempfile
import time

import skbio
from q2_types.feature_data import BLAST6Format, DNAFASTAFormat

from ._alignment_evidence import (
    HIT_COLUMNS, QUERY_COLUMNS, inventory_fingerprint, validate_hit)
from ._cluster_features import run_command
from ._format import MetaVsearchReferenceAlignmentsDirFmt


USER_FIELDS = (
    'query', 'target', 'qstrand', 'ids', 'alnlen', 'mism', 'gaps', 'opens',
    'exts', 'ql', 'tl', 'qilo', 'qihi', 'tilo', 'tihi', 'raw', 'caln')


def _stage_fasta(source, destination, conn, kind):
    count = 0
    maximum = 1
    with open(destination, 'w', encoding='utf-8') as fh:
        for count, seq in enumerate(skbio.io.read(
                str(source), format='fasta', constructor=skbio.DNA), 1):
            original = seq.metadata['id']
            sequence = str(seq).upper()
            if (not sequence or not original or
                    any(char in original for char in '\t\r\n')):
                raise ValueError('FASTA IDs and sequences must be nonempty '
                                 'and IDs must not contain tabs/newlines.')
            token = '%s%012d' % (kind, count)
            digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
            try:
                conn.execute('INSERT INTO records VALUES '
                             '(?, ?, ?, ?, ?, ?, 0, 0, 1)',
                             (token, kind, original, len(sequence), digest,
                              count))
            except sqlite3.IntegrityError as exc:
                raise ValueError('Duplicate %s sequence identifier: %s' %
                                 (kind, original)) from exc
            fh.write('>%s\n%s\n' % (token, sequence))
            maximum = max(maximum, len(sequence))
    if not count:
        raise ValueError('Input sequence collections must not be empty.')
    conn.commit()

    def fingerprint(order):
        return inventory_fingerprint(conn.execute(
            'SELECT original_id, length, sha256 FROM records '
            'WHERE kind=? ORDER BY ' + order, (kind,)))
    return {'count': count, 'sha256': fingerprint('original_id'),
            'order_sha256': fingerprint('ordinal'), 'max_length': maximum}


def _search_parameters(identity_floor, strand, search_mode, threads, maximum):
    if not 0 < identity_floor <= 1:
        raise ValueError('identity_floor must be >0 and at most 1.')
    if strand not in ('plus', 'both'):
        raise ValueError('strand must be plus or both.')
    if search_mode not in ('exhaustive', 'heuristic'):
        raise ValueError('search_mode must be exhaustive or heuristic.')
    if not isinstance(threads, int) or not 0 <= threads <= 1024:
        raise ValueError('threads must be an integer from 0 to 1024.')
    exhaustive = search_mode == 'exhaustive'
    return {'id': identity_floor, 'strand': strand, 'threads': threads,
            'iddef': 2, 'qmask': 'none', 'dbmask': 'dust',
            'match': 2, 'mismatch': -4, 'gapopen': '20I/2E',
            'gapext': '2I/1E', 'wordlength': 8, 'minseqlength': 1,
            'maxseqlength': maximum, 'query_cov': 0.0, 'target_cov': 0.0,
            'maxaccepts': 0 if exhaustive else 1,
            'maxrejects': 0 if exhaustive else 32,
            'minwordmatches': 0 if exhaustive else 12, 'maxhits': 0,
            'fasta_width': 0}


def _parse_hits(path, conn):
    with open(path, encoding='utf-8') as fh:
        for line_number, line in enumerate(fh, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) != len(USER_FIELDS):
                raise ValueError('Invalid VSEARCH userout at line %d.' %
                                 line_number)
            query = conn.execute('SELECT original_id, length, best_matches, '
                                 'best_length FROM records WHERE token=? '
                                 'AND kind="q"', (fields[0],)).fetchone()
            reference = conn.execute('SELECT original_id, length FROM '
                                     'records WHERE token=? AND kind="r"',
                                     (fields[1],)).fetchone()
            if query is None or reference is None:
                raise ValueError('VSEARCH returned an unknown sequence ID.')
            values = list(map(int, fields[3:16]))
            matches, length, mismatches = values[:3]
            if values[6:8] != [query[1], reference[1]]:
                raise ValueError('VSEARCH returned inconsistent lengths.')
            row = dict(zip(HIT_COLUMNS[:17],
                           [query[0], reference[0], fields[2]] + values +
                           [fields[16]]))
            row.update(identity=matches / length,
                       query_coverage=(matches + mismatches) / query[1],
                       reference_coverage=(
                           matches + mismatches) / reference[1],
                       is_exact_best='0')
            validate_hit(row)
            try:
                conn.execute('INSERT INTO hits VALUES (' +
                             ','.join('?' for _ in HIT_COLUMNS[:17]) + ')',
                             tuple(row[key] for key in HIT_COLUMNS[:17]))
            except sqlite3.IntegrityError as exc:
                raise ValueError('Duplicate query/reference/strand hit.') \
                    from exc
            conn.execute('UPDATE records SET hit_count=hit_count+1 '
                         'WHERE token=?', (fields[0],))
            if matches * query[3] > query[2] * length:
                conn.execute('UPDATE records SET best_matches=?, '
                             'best_length=? WHERE token=?',
                             (matches, length, fields[0]))
    conn.commit()


def _write_outputs(conn, evidence, blast6):
    root = Path(str(evidence))
    with open(root / 'queries.tsv', 'w', encoding='utf-8', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(QUERY_COLUMNS)
        for row in conn.execute('SELECT original_id, length, sha256, '
                                'hit_count FROM records WHERE kind="q" '
                                'ORDER BY original_id'):
            writer.writerow((*row, 'matched' if row[3] else 'no_hit'))
    with open(root / 'hits.tsv', 'w', encoding='utf-8', newline='') as fh, \
            open(str(blast6), 'w', encoding='utf-8', newline='') as bf:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        blast_writer = csv.writer(bf, delimiter='\t', lineterminator='\n')
        writer.writerow(HIT_COLUMNS)
        for values in conn.execute('SELECT * FROM hits ORDER BY '
                                   'query_id, reference_id, strand'):
            row = dict(zip(HIT_COLUMNS[:17], values))
            matches, length = row['matches'], row['alignment_length']
            best = conn.execute('SELECT best_matches, best_length FROM '
                                'records WHERE kind="q" AND original_id=?',
                                (row['query_id'],)).fetchone()
            pairs = matches + row['mismatches']
            derived = [format(matches / length, '.17g'),
                       format(pairs / row['query_length'], '.17g'),
                       format(pairs / row['reference_length'], '.17g'),
                       int(matches * best[1] == best[0] * length)]
            writer.writerow((*values, *derived))
            reverse = row['strand'] == '-'
            blast_writer.writerow([
                row['query_id'], row['reference_id'],
                format(100 * matches / length, '.1f'), length,
                row['mismatches'], row['gap_opens'],
                row['query_length'] if reverse else 1,
                1 if reverse else row['query_length'],
                1, row['reference_length'], -1, 0])
        # Explicit no-hit records keep standard BLAST6 valid even when all
        # queries fail. These are never included in the evidence hit table.
        for query, in conn.execute('SELECT original_id FROM records '
                                   'WHERE kind="q" AND hit_count=0 '
                                   'ORDER BY original_id'):
            blast_writer.writerow([query, '*', '0.0', 0, 0, 0, 0, 0, 0, 0,
                                   -1, 0])


def align_closed_reference(sequences: DNAFASTAFormat,
                           reference_sequences: DNAFASTAFormat,
                           identity_floor: float = 0.90,
                           strand: str = 'plus',
                           search_mode: str = 'exhaustive',
                           threads: int = 1
                           ) -> (MetaVsearchReferenceAlignmentsDirFmt,
                                 BLAST6Format):
    evidence = MetaVsearchReferenceAlignmentsDirFmt()
    blast6 = BLAST6Format()
    with tempfile.TemporaryDirectory(prefix='meta-vsearch-') as tmp:
        root = Path(tmp)
        with closing(sqlite3.connect(root / 'alignment.sqlite')) as conn:
            conn.execute('CREATE TABLE records (token TEXT PRIMARY KEY, '
                         'kind TEXT, original_id TEXT, length INTEGER, '
                         'sha256 TEXT, ordinal INTEGER, hit_count INTEGER, '
                         'best_matches INTEGER, best_length INTEGER, '
                         'UNIQUE(kind, original_id))')
            columns = ', '.join(key + (' TEXT' if key in (
                'query_id', 'reference_id', 'strand', 'cigar') else ' INTEGER')
                for key in HIT_COLUMNS[:17])
            conn.execute('CREATE TABLE hits (' + columns +
                         ', PRIMARY KEY(query_id, reference_id, strand))')
            queries = _stage_fasta(
                sequences, root / 'queries.fasta', conn, 'q')
            references = _stage_fasta(reference_sequences,
                                      root / 'references.fasta', conn, 'r')
            params = _search_parameters(
                identity_floor, strand, search_mode, threads,
                max(queries['max_length'], references['max_length']))
            cmd = ['vsearch', '--usearch_global', str(root / 'queries.fasta'),
                   '--db', str(root / 'references.fasta'),
                   '--userout', str(root / 'hits.userout'),
                   '--userfields', '+'.join(USER_FIELDS),
                   '--log', str(root / 'vsearch.log')]
            for key, value in params.items():
                cmd.extend(['--' + key, str(value)])
            started = time.perf_counter()
            try:
                run_command(cmd)
            except subprocess.CalledProcessError as exc:
                raise RuntimeError('VSEARCH alignment failed with exit code '
                                   '%s.' % exc.returncode) from exc
            elapsed = time.perf_counter() - started
            log = (root / 'vsearch.log').read_text(encoding='utf-8')
            version = re.search(r'^vsearch v(\d+)\.(\d+)\.(\d+)[^\n]*',
                                log, re.MULTILINE)
            if (version is None or
                    tuple(map(int, version.groups())) < (2, 3, 1)):
                raise ValueError('VSEARCH >=2.3.1 with a versioned log is '
                                 'required for a reliable zero-word bypass.')
            _parse_hits(root / 'hits.userout', conn)
            _write_outputs(conn, evidence, blast6)
            manifest = {
                'schema_version': 1, 'vsearch_version': version.group(0),
                'search_mode': search_mode,
                'search_complete': search_mode == 'exhaustive',
                'completeness_scope': 'all input references and configured '
                                      'query strands; one alignment per pair '
                                      'and strand under VSEARCH scoring',
                'identity_definition': 2, 'identity_floor': identity_floor,
                'exact_best_definition': 'maximum matches/alignment_length; '
                                         'integer cross multiplication',
                'identity_tolerance': 0,
                'strand': strand, 'effective_parameters': params,
                'effective_flags': {'top_hits_only': False, 'hardmask': False,
                                    'n_mismatch': False, 'sizein': False,
                                    'self': False, 'selfid': False},
                'queries': queries, 'references': references,
                'fingerprint_algorithm': 'sha256-json-records-v1; uppercase '
                                         'sequence; sorted original IDs',
                'retained_hits': 'all accepted alignments',
                'vsearch_elapsed_seconds': elapsed,
                'vsearch_log': log.replace(str(root), '<workdir>'),
                'command': [part.replace(str(root), '<workdir>')
                            for part in cmd]
            }
            (Path(str(evidence)) / 'search.json').write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + '\n',
                encoding='utf-8')
    return evidence, blast6
