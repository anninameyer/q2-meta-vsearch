# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import csv
from pathlib import Path

import pandas as pd
import qiime2
from q2_types.feature_data import (
    SequenceCharacteristicsDirectoryFormat, TSVTaxonomyDirectoryFormat)

from .plugin_setup import plugin
from ._format import MetaVsearchUchimeStatsFmt

# many of the numeric fields will contain * if a query is
# not chimeric, so being safe and making all fields other than
# score strings
_uchime_stats_header = collections.OrderedDict([
     ('score', float),
     ('feature-id', str),
     ('A', str),
     ('B', str),
     ('T', str),
     ('idQM', str),
     ('idQA', str),
     ('idQB', str),
     ('idAB', str),
     ('idQT', str),
     ('LY', str),
     ('LN', str),
     ('LA', str),
     ('RY', str),
     ('RN', str),
     ('RA', str),
     ('div', str),
     ('YN', str)])


def _stats_to_df(ff):
    df = pd.read_csv(str(ff), sep='\t', index_col='feature-id',
                     names=_uchime_stats_header.keys(),
                     dtype=_uchime_stats_header)
    return df


@plugin.register_transformer
def _1(ff: MetaVsearchUchimeStatsFmt) -> qiime2.Metadata:
    return qiime2.Metadata(_stats_to_df(ff))


@plugin.register_transformer
def _2(ff: MetaVsearchUchimeStatsFmt) -> pd.DataFrame:
    return _stats_to_df(ff)


def _metadata_directory_view(path):
    """Read IDs as strings; honor metadata directives when present.

    These directory-level routes supplement, rather than replace, q2-types'
    file-level transformers, which infer numeric IDs and pandas NA spellings.
    They apply to the existing canonical types without registering new formats.
    """
    with open(path, encoding='utf-8', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        header = next(reader)
        following = next(reader, [])
    if following and following[0].startswith('#q2:'):
        return qiime2.Metadata.load(str(path))
    frame = pd.read_csv(path, sep='\t', dtype={header[0]: str},
                        keep_default_na=False,
                        na_values={column: [''] for column in header[1:]})
    frame = frame.set_index(header[0])
    frame.index.name = 'id'
    for column in frame.select_dtypes(include='bool'):
        frame[column] = frame[column].astype(str)
    return qiime2.Metadata(frame)


@plugin.register_transformer
def _3(ff: SequenceCharacteristicsDirectoryFormat) -> qiime2.Metadata:
    path = Path(str(ff)) / 'sequence_characteristics.tsv'
    return _metadata_directory_view(path)


@plugin.register_transformer
def _4(ff: SequenceCharacteristicsDirectoryFormat) -> pd.DataFrame:
    return _3(ff).to_dataframe()


@plugin.register_transformer
def _5(ff: TSVTaxonomyDirectoryFormat) -> qiime2.Metadata:
    return _metadata_directory_view(Path(str(ff)) / 'taxonomy.tsv')


@plugin.register_transformer
def _6(ff: TSVTaxonomyDirectoryFormat) -> pd.DataFrame:
    frame = _5(ff).to_dataframe()
    frame.index.name = 'Feature ID'
    return frame


@plugin.register_transformer
def _7(ff: TSVTaxonomyDirectoryFormat) -> pd.Series:
    return _6(ff)['Taxon']
