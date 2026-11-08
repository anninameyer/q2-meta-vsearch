# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

"""QIIME 2 assignment adapter; never executes VSEARCH."""

import csv

import biom
import pandas as pd
import qiime2
from q2_types.feature_data import (
    SequenceCharacteristicsFormat, TSVTaxonomyFormat)

from ._assignment import (
    assign_reference_alignments, NUMERIC_COLUMNS, SUMMARY_COLUMNS)
from ._format import MetaVsearchReferenceAlignmentsDirFmt


def read_class_weights(table):
    if table.shape[1] != 1:
        raise ValueError('class_weights must contain exactly one weight '
                         'vector (sample); its ID need not be Weight.')
    return dict(zip(table.ids(axis='observation'),
                    table.data(table.ids()[0], axis='sample')))


def read_reference_taxonomy(path):
    taxonomy = {}
    with open(path, encoding='utf-8', newline='') as fh:
        rows = (row for row in csv.reader(fh, delimiter='\t') if row)
        header = next(rows, [])
        if header[:2] != ['Feature ID', 'Taxon']:
            raise ValueError('Reference taxonomy needs Feature ID and Taxon.')
        for fields in rows:
            if len(fields) != len(header):
                raise ValueError('Malformed reference taxonomy row.')
            reference = fields[0]
            if not reference or reference in taxonomy:
                raise ValueError('Empty or duplicate reference taxonomy ID.')
            taxonomy[reference] = fields[1]
    return taxonomy


def assign_closed_reference(alignments: MetaVsearchReferenceAlignmentsDirFmt,
                            reference_taxonomy: TSVTaxonomyFormat,
                            identity_threshold: float = 0.90,
                            identity_tolerance: float = 0.0,
                            unresolved_policy: str = 'deterministic-reference',
                            max_explicit_species: int = 2,
                            class_weights: biom.Table = None,
                            min_weight_difference: float = 0.0,
                            min_weight_ratio: float = 1.0,
                            missing_weight_policy: str = 'unresolved'
                            ) -> (SequenceCharacteristicsFormat,
                                  TSVTaxonomyFormat):
    alignments.validate('max')
    result = assign_reference_alignments(
        str(alignments), read_reference_taxonomy(str(reference_taxonomy)),
        identity_threshold, identity_tolerance, unresolved_policy,
        max_explicit_species,
        None if class_weights is None else read_class_weights(class_weights),
        min_weight_difference, min_weight_ratio, missing_weight_policy)
    frame = pd.DataFrame(result.summary).set_index('id')
    frame = frame.loc[:, SUMMARY_COLUMNS]
    for column in NUMERIC_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors='raise')
    metadata = qiime2.Metadata(frame)
    summary = SequenceCharacteristicsFormat()
    # Type directives keep numeric-looking reference IDs categorical and
    # distinguish absent statistics from categorical strings such as "NA".
    metadata.save(str(summary))
    taxonomy = TSVTaxonomyFormat()
    with open(str(taxonomy), 'w', encoding='utf-8', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
        writer.writerow(['Feature ID', 'Taxon'])
        writer.writerows(result.taxonomy.items())
    return summary, taxonomy
