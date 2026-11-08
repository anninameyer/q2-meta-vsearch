# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

"""Conservative lineage consensus. Reference multiplicity supplies no votes."""

from dataclasses import dataclass
import re


RANKS = ('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
PREFIXES = ('d', 'p', 'c', 'o', 'f', 'g', 's')
DOMAINS = {'Bacteria', 'Archaea', 'Eukaryota', 'Eukarya', 'Viruses'}
UNKNOWN = {'', 'unassigned', 'unclassified', 'uncultured', 'unknown',
           'unidentified', 'metagenome', 'environmental sample', 'na', 'n/a',
           'none', 'null', 'nan', '-', '?'}


@dataclass(frozen=True)
class Lineage:
    # A contiguous, explicitly supported prefix of the seven-rank hierarchy.
    values: tuple
    root_prefix: str = 'd'

    @property
    def taxonomy(self):
        prefixes = (self.root_prefix, *PREFIXES[1:])
        return '; '.join('%s__%s' % (prefix, value)
                         for prefix, value in zip(prefixes, self.values)) \
            or 'Unassigned'


def _known(value):
    lower = value.casefold()
    return (lower not in UNKNOWN and
            not lower.startswith(('uncultured ', 'unclassified ',
                                  'unidentified ', 'unknown ')))


def _species(value, genus):
    # Require an explicit genus and a conservative binomial/epithet spelling.
    # Do not turn sp., strain labels, or arbitrary organism names into species.
    words = value.split()
    if len(words) == 1:
        epithet = words[0]
    elif len(words) == 2 and words[0] == genus:
        epithet = words[1]
    else:
        return None
    if (not re.fullmatch(r'[A-Z][A-Za-z_-]+', genus) or
            not re.fullmatch(r'[a-z][a-z-]+', epithet) or
            not _known(epithet) or
            epithet in {'sp', 'spp', 'bacterium', 'archaeon', 'organism'}):
        return None
    return genus + ' ' + epithet


def parse_lineage(taxonomy):
    if taxonomy is None or not isinstance(taxonomy, str):
        return Lineage(())
    tokens = [part.strip() for part in taxonomy.split(';')]
    # Allow a conventional trailing semicolon, without filling missing ranks.
    while tokens and not tokens[-1]:
        tokens.pop()
    values = []
    root_prefix = 'd'
    for index, token in enumerate(tokens[:7]):
        prefixed = re.fullmatch(r'([dkpcofgs])__(.*)', token)
        numbered = re.fullmatch(r'D_([0-6])__(.*)', token)
        if prefixed:
            prefix, value = prefixed.groups()
            rank = 0 if prefix == 'k' else PREFIXES.index(prefix)
            if (prefix == 'k' and index == 0 and
                    value.strip() not in DOMAINS):
                root_prefix = 'k'
        elif numbered:
            rank, value = int(numbered[1]), numbered[2]
        elif '__' not in token and tokens[0] in DOMAINS:
            rank, value = index, token
        else:
            break
        value = value.strip()
        if rank != index or not _known(value):
            break
        if rank == 6:
            value = _species(value, values[5])
            if value is None:
                break
        values.append(value)
    return Lineage(tuple(values), root_prefix)


def consensus_taxonomy(lineages, missing_count, max_explicit_species):
    """Return LCA, readable label, rank, status, distinct known lineage count.

    Missing best-hit annotations veto a claimed common ancestor: an unknown
    candidate could belong to any lineage. Identical incomplete annotations
    support their stated common prefix, without asserting missing ranks.
    """
    distinct = {lineage for lineage in lineages if lineage.values}
    if missing_count:
        status = 'partial_missing_taxonomy' if distinct else 'missing_taxonomy'
        return 'Unassigned', 'Unassigned', 'Unassigned', status, len(distinct)
    shared = []
    roots = {lineage.root_prefix for lineage in distinct}
    for values in zip(*(lineage.values for lineage in distinct)):
        if len(roots) != 1:
            break
        if len(set(values)) != 1:
            break
        shared.append(values[0])
    root = next(iter(roots)) if len(roots) == 1 else 'd'
    lca = Lineage(tuple(shared), root).taxonomy
    rank = ('Kingdom' if len(shared) == 1 and root == 'k' else
            RANKS[len(shared) - 1] if shared else 'Unassigned')
    status = ('ambiguous' if len(distinct) > 1 else
              'resolved' if len(shared) == 7 else 'incomplete')
    label = shared[-1] if shared else 'Unassigned'
    # Compact labels require complete, validated species in the same lineage.
    if len(shared) == 6 and all(len(item.values) == 7 for item in distinct):
        genus = shared[5]
        epithets = sorted({item.values[6].split(' ', 1)[1]
                           for item in distinct})
        label = (genus + ' ' + '|'.join(epithets)
                 if len(epithets) <= max_explicit_species else genus + ' ssg.')
    return lca, label, rank, status, len(distinct)
