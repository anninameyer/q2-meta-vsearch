# ----------------------------------------------------------------------------
# Copyright (c) 2026, q2-meta-vsearch developers.
# Distributed under the terms of the Modified BSD License.
# ----------------------------------------------------------------------------

"""Taxon-level priors applied exclusively to exact sequence ties."""

from fractions import Fraction
import math

from ._alignment_evidence import identity_fraction
from ._taxonomy import consensus_taxonomy, parse_lineage


WEIGHT_COLUMNS = (
    'reference_with_weights', 'taxonomy_with_weights', 'selected_weight',
    'second_best_weight', 'weight_difference', 'weight_ratio',
    'n_weighted_best_taxa', 'n_weighted_best_reference_hits',
    'weight_resolved', 'assignment_changed', 'taxonomy_changed',
    'missing_weight_count',
    'weight_assignment_status')
WEIGHT_NUMERIC_COLUMNS = {
    'selected_weight', 'second_best_weight', 'weight_difference',
    'weight_ratio',
    'n_weighted_best_taxa', 'n_weighted_best_reference_hits',
    'missing_weight_count'}


def normalize_taxonomy_key(value):
    """Whitespace only: never change rank prefixes, case, or missing ranks."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError('Weight taxonomy keys must be nonempty strings.')
    return ';'.join(' '.join(part.split()) for part in value.split(';'))


def validate_weight_parameters(min_weight_difference, min_weight_ratio,
                               missing_weight_policy):
    for name, value, lower in (
            ('min_weight_difference', min_weight_difference, 0),
            ('min_weight_ratio', min_weight_ratio, 1)):
        if (isinstance(value, bool) or not math.isfinite(value) or
                value < lower):
            raise ValueError('%s must be finite and >= %s.' % (name, lower))
    if missing_weight_policy not in ('unresolved', 'error'):
        raise ValueError('Invalid missing_weight_policy.')


def prepare_weights(class_weights):
    """Validate an ordinary mapping, rejecting even equal-valued aliases."""
    normalized = {}
    originals = {}
    for key, value in class_weights.items():
        canonical = normalize_taxonomy_key(key)
        if canonical in normalized:
            raise ValueError('Taxonomy normalization collision between %r '
                             'and %r.' % (originals[canonical], key))
        if (isinstance(value, bool) or not math.isfinite(value) or value < 0):
            raise ValueError('Ecological weight for %r must be finite and '
                             'nonnegative.' % key)
        normalized[canonical] = Fraction(str(value))
        originals[canonical] = key
    if not normalized:
        raise ValueError('class_weights must contain at least one taxon.')
    return normalized


def without_weights(row):
    """Expose equivalent decisions and absent measurements without priors."""
    row.update(dict.fromkeys(WEIGHT_COLUMNS))
    row.update(reference_with_weights=row['reference_without_weights'],
               taxonomy_with_weights=row['taxonomy_without_weights'],
               n_weighted_best_taxa=row['n_best_taxa'],
               n_weighted_best_reference_hits=row['n_best_reference_hits'],
               weight_resolved='false', assignment_changed='false',
               taxonomy_changed='false', missing_weight_count=0,
               weight_assignment_status='not_provided')


def weight_query(fetch_hits, taxonomy, row, weights,
                 min_weight_difference, min_weight_ratio,
                 missing_weight_policy, unresolved_policy,
                 max_explicit_species):
    """Extend a baseline row; return final taxonomy and selected reference.

    No sum over reference sequences and no candidate-wise renormalization.
    Integer-ratio equality is rechecked before accessing any weight.
    """
    row.update(dict.fromkeys(WEIGHT_COLUMNS))
    row.update(n_weighted_best_taxa=0, n_weighted_best_reference_hits=0,
               missing_weight_count=0, weight_resolved='false',
               assignment_changed='false', taxonomy_changed='false',
               taxonomy_with_weights=row['taxonomy_without_weights'],
               weight_assignment_status='no_match')
    if row['matched'] == 'false':
        return row['taxonomy_with_weights'], None

    maximum = Fraction(row['best_identity_numerator'],
                       row['best_identity_denominator'])
    candidates, unknown = {}, set()
    for hit in fetch_hits():
        if identity_fraction(hit) != maximum:
            continue
        ref = hit['reference_id']
        raw = taxonomy.get(ref)
        if not parse_lineage(raw).values:
            unknown.add(ref)
            continue
        key = normalize_taxonomy_key(raw)
        candidates.setdefault(key, set()).add(ref)
    missing = sorted(key for key in candidates if key not in weights)
    row['missing_weight_count'] = len(missing) + len(unknown)
    remaining = set(candidates)
    status = 'already_resolved'
    if missing or unknown:
        if missing_weight_policy == 'error':
            raise ValueError('Missing ecological weights for candidate '
                             'lineages %r; missing/unusable reference '
                             'taxonomy for IDs %r.' %
                             (missing, sorted(unknown)))
        status = 'missing_taxonomy' if unknown else 'missing_weights'
    else:
        ordered = sorted(candidates, key=lambda key: (-weights[key], key))
        top = weights[ordered[0]]
        if len(ordered) > 1:
            second = weights[ordered[1]]
            difference = top - second
            row.update(second_best_weight=float(second),
                       weight_difference=float(difference))
            # Positive/zero is mathematically infinite. Leave the numeric
            # diagnostic missing rather than serialize infinity as metadata.
            ratio = (top / second if second else
                     Fraction(1) if not top else None)
            if ratio is not None:
                try:
                    numeric_ratio = float(ratio)
                except OverflowError:
                    numeric_ratio = math.inf
                row['weight_ratio'] = (numeric_ratio if
                                       math.isfinite(numeric_ratio) else None)
            decisive = (
                top > second and
                difference >= Fraction(str(min_weight_difference)) and
                top >= second * Fraction(str(min_weight_ratio)))
            if decisive:
                remaining = {ordered[0]}
                status = 'resolved'
            else:
                status = ('equal_weights' if top == second else
                          'insufficient_margin')

    refs = unknown | set().union(*(candidates[key] for key in remaining))
    lca, _, _, tax_status, _ = consensus_taxonomy(
        {parse_lineage(taxonomy[ref]) for ref in refs - unknown}, len(unknown),
        max_explicit_species)
    excluded = (unresolved_policy == 'exclude' and tax_status in (
        'ambiguous', 'missing_taxonomy', 'partial_missing_taxonomy'))
    selected = None if excluded else min(refs)
    if len(remaining) == 1 and not unknown and not missing:
        row['selected_weight'] = float(weights[next(iter(remaining))])
    row.update(
        reference_with_weights=selected, taxonomy_with_weights=lca,
        n_weighted_best_taxa=len(remaining),
        n_weighted_best_reference_hits=len(refs),
        weight_assignment_status=status,
        weight_resolved=str(
            status == 'resolved' and
            row['taxonomic_assignment_status'] == 'ambiguous' and
            tax_status != 'ambiguous').lower(),
        assignment_changed=str(
            selected != row['reference_without_weights']).lower(),
        taxonomy_changed=str(lca != row['taxonomy_without_weights']).lower())
    return lca, selected
