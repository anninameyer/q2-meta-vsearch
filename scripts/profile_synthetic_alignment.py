#!/usr/bin/env python
"""Profile one committed software-validation fixture in a fresh process.

These synthetic cases measure software behavior only, not microbiome accuracy.
Run once per case so child-process memory high-water marks are independent.
"""

import argparse
import csv
import json
from pathlib import Path
import resource
import shutil
import sys
import time

from q2_types.feature_data import DNAFASTAFormat
from q2_meta_vsearch._align import align_closed_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['basic', 'precision', 'geometry',
                                           'suboptimal'], default='basic')
    parser.add_argument('--search-mode', choices=['exhaustive', 'heuristic'],
                        default='exhaustive')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    fixtures = Path(__file__).resolve().parents[1] / (
        'q2_meta_vsearch/tests/data/alignment')
    queries, references = {
        'basic': ('queries.fasta', 'references.fasta'),
        'precision': ('precision-query.fasta', 'precision-references.fasta'),
        'geometry': ('geometry.fasta', 'reference-single.fasta'),
        'suboptimal': ('suboptimal-query.fasta', 'suboptimal-references.fasta')
    }[args.case]
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    evidence, blast6 = align_closed_reference(
        DNAFASTAFormat(fixtures / queries, mode='r'),
        DNAFASTAFormat(fixtures / references, mode='r'),
        strand='both' if args.case == 'geometry' else 'plus',
        search_mode=args.search_mode)
    evidence.validate('max')
    elapsed = time.perf_counter() - started
    manifest = json.loads((Path(str(evidence)) / 'search.json').read_text())
    with open(Path(str(evidence)) / 'hits.tsv') as fh:
        count = sum(1 for _ in csv.DictReader(fh, delimiter='\t'))
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    report = {
        'synthetic_software_validation_only': True,
        'case': args.case, 'search_mode': args.search_mode,
        'queries': manifest['queries']['count'],
        'references': manifest['references']['count'],
        'accepted_alignments': count,
        'vsearch_version': manifest['vsearch_version'],
        'vsearch_elapsed_seconds': manifest['vsearch_elapsed_seconds'],
        'alignment_and_validation_seconds': elapsed,
        'child_process_peak_rss_bytes': usage.ru_maxrss * (
            1 if sys.platform == 'darwin' else 1024),
        'memory_scope': 'fresh-process RUSAGE_CHILDREN high-water mark; '
                        'not parent Python memory or a scaling benchmark',
        'query_fingerprint': manifest['queries']['sha256'],
        'reference_fingerprint': manifest['references']['sha256']}
    shutil.copytree(str(evidence), args.output / 'evidence')
    shutil.copyfile(str(blast6), args.output / 'blast6.tsv')
    (args.output / 'profile.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
