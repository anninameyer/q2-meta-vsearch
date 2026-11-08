# Ecological prior weighting (M5)

`assign-closed-reference` optionally consumes a standard
`FeatureTable[RelativeFrequency]` artifact, including the output of
`qiime clawback generate-class-weights`. q2-clawback is not a runtime dependency.
There is no OTU aggregation, new reassignment action, or combined pipeline in M5.

## Action signature

```text
assign_closed_reference(
    alignments: FeatureData[ReferenceAlignments],
    reference_taxonomy: FeatureData[Taxonomy],
    identity_threshold: Float = 0.90,
    identity_tolerance: Float = 0.0,
    unresolved_policy: Choices['deterministic-reference', 'exclude']
        = 'deterministic-reference',
    max_explicit_species: Int = 2,
    class_weights: FeatureTable[RelativeFrequency] = None,
    min_weight_difference: Float = 0.0,
    min_weight_ratio: Float = 1.0,
    missing_weight_policy: Choices['unresolved', 'error'] = 'unresolved'
) -> (
    assignment_metadata: FeatureData[SequenceCharacteristics],
    taxonomy: FeatureData[Taxonomy]
)
```

```bash
qiime meta-vsearch assign-closed-reference \
  --i-alignments alignments.qza \
  --i-reference-taxonomy reference-taxonomy.qza \
  --i-class-weights clawback-weights.qza \
  --p-identity-threshold 0.97 \
  --p-min-weight-difference 0.0 \
  --p-min-weight-ratio 1.0 \
  --p-missing-weight-policy unresolved \
  --o-assignment-metadata weighted-assignments.qza \
  --o-taxonomy weighted-taxonomy.qza
```

Omit `--i-class-weights` for the original unweighted decisions. No ASV abundance
table is required. The action never invokes VSEARCH or modifies the evidence.
Passing another vector reuses the same evidence within this action; the complete
M7 reassignment workflow remains future work.

## Input validation and lineage matching

The weight table must have exactly **one sample/vector**. Its sample identifier
may be `Weight` or any other identifier. Features are complete taxonomy strings;
values must be finite and nonnegative. Explicit zeros are valid; an absent
feature is missing, never zero. Empty vectors, multiple samples, negative
values, NaN, and infinity are rejected. Weights are used as supplied without
per-ASV renormalization. This plugin does not repair or rescale a vector.

Mapping uses the complete supplied reference annotation, including incomplete
or placeholder ranks. It does not use only the species name or the shorter
parsed prefix used for conservative LCA reporting. Normalization strips leading
and trailing whitespace on each semicolon-delimited component, collapses runs
of whitespace within components, and ignores spacing around semicolons. It
preserves case, rank prefixes, empty components, and trailing semicolons.
For example, `d__Bacteria; p__Bacillota` matches
` d__Bacteria ;   p__Bacillota `.

Distinct weight feature IDs that normalize to the same key are a **collision**
and raise an error even if their numeric weights agree. Multiple references
sharing that key are expected and receive the same single taxon weight.
Different rank spellings (`d__` versus `k__`), synonyms, database versions,
terminal-species-only labels, or omitted ranks are not reconciled. A mismatch
is reported through missing weights; `error` makes it fatal. Database identities
cannot be verified from taxonomy labels alone: users must supply compatible
reference taxonomy and ecological priors.

A best reference without a usable taxonomy root cannot be mapped safely. It is
reported as missing taxonomy and blocks weighting. Identical partial annotations
may be weighted as supplied, but never gain an inferred rank or species name.
The conservative M3–M4 taxonomy parser and LCA rules remain unchanged.

## Exact decision and margin rules

1. Compute the maximum sequence identity as an integer-count ratio, using all
   stored accepted alignments. Apply the assignment threshold as before.
2. Preserve the original unweighted decision and all original metadata.
3. Consider **only exact-best reference sequences**. Lower-identity and near-best
   diagnostic hits cannot supply candidates, taxonomy, or weights.
4. Group the candidates by their normalized complete input taxonomy. Look up
   one weight per group, without summing over references or multiplying by their
   count. Ten references from A supply the same support as one reference from A.
5. With complete weights, compare the largest taxon weight `w1` and runner-up
   `w2`. A winner requires all three conditions:

   ```text
   w1 > w2
   w1 - w2 >= min_weight_difference
   w1 >= w2 * min_weight_ratio
   ```

   Both margins are inclusive; the top itself must be unique. Defaults 0 and 1
   allow any unique maximum. Parameters must be finite, difference >= 0 and
   ratio >= 1. Decimal representations of supplied numbers are compared as exact
   fractions, avoiding subtraction rounding at a configured margin boundary.
6. If decisive, retain every reference in the winning taxon. If the top is tied
   or either margin fails, conservatively retain **all original exact-best
   candidate taxa**, including lower-weight taxa. There is no partial pruning.
   A single candidate taxon needs no contest and is `already_resolved`.
7. Compute the hierarchical LCA of the remaining references. Select the
   lexicographically first original reference ID within that set. Under
   `exclude`, suppress selection if the final taxonomy remains ambiguous or
   missing. A weighted winner may therefore restore a reference excluded by
   the unweighted decision.

For positive `w1` and zero `w2`, the mathematical ratio is infinite and every
finite ratio requirement passes; the absolute difference must still pass.
Zero/zero is an unresolved tie with diagnostic ratio 1. Ratios use exact
multiplication for decisions, including extremely small positive weights.
The numeric `weight_ratio` is missing for positive/zero or an overflow beyond
finite floating-point range. Its absence is not evidence of a failed criterion.

When any exact-best candidate lacks a weight, `missing_weight_policy=unresolved`
retains the entire candidate set without ranking the available subset. Missing
reference taxonomy behaves the same way. Reference selection follows the
original deterministic/exclude policy and taxonomic ambiguity remains explicit.
The `error` policy raises with the affected lineages and reference IDs. Missing
weights on lower-identity hits do not matter. A missing weight does not make an
otherwise unambiguous single annotated taxon ambiguous, but it is still reported.

## Weighted metadata schema

The original 30 columns and their values remain unweighted; see the
[original schema](assignment.md). In particular, `n_best_taxa`, `taxonomy_lca`,
`taxonomy_ambiguity_label`, `taxonomic_assignment_status`, and
`taxonomic_resolution_rank` describe the original exact-best consensus.
Sequence statistics still describe `statistics_reference`, even if a prior
selects a different reference with different coverage. The full alignment
artifact retains those alternative statistics.

These 13 columns are added, making 43 columns in both weighted and unweighted
outputs. All original ASV identifiers appear exactly once, including unmatched
ASVs. Boolean fields are categorical `true`/`false`; blank numeric cells are
missing values, never imputed zero weights.

| Column | Metadata type | Meaning with priors |
| --- | --- | --- |
| `reference_with_weights` | categorical | Final deterministic reference, or blank if excluded/unmatched |
| `taxonomy_with_weights` | categorical | Final hierarchical LCA, identical to the primary taxonomy artifact |
| `selected_weight` | numeric | Weight of the sole remaining taxon when all candidate weights/taxonomy are present; blank for unresolved contests or missing evidence |
| `second_best_weight` | numeric | Runner-up among original distinct exact-best taxon weights; blank for a single taxon or incomplete weights |
| `weight_difference` | numeric | Original top minus runner-up; reported even if margins fail |
| `weight_ratio` | numeric | Top/runner-up; 1 for zero/zero, blank for positive/zero or overflow; blank without a complete contest |
| `n_weighted_best_taxa` | integer-valued numeric | Remaining normalized complete taxonomy keys; unknown annotations are not invented taxa |
| `n_weighted_best_reference_hits` | integer-valued numeric | Remaining exact-best reference IDs, even if selection is excluded |
| `weight_resolved` | categorical boolean | True when a decisive prior removes original taxonomic disagreement; a reference change alone is insufficient |
| `assignment_changed` | categorical boolean | Final reference differs from `reference_without_weights`, including a change to/from exclusion |
| `taxonomy_changed` | categorical boolean | Final hierarchical taxonomy string differs from `taxonomy_without_weights` |
| `missing_weight_count` | integer-valued numeric | Number of distinct known candidate keys absent from the vector, plus the number of best reference IDs without usable taxonomy |
| `weight_assignment_status` | categorical | `resolved`, `already_resolved`, `equal_weights`, `insufficient_margin`, `missing_weights`, `missing_taxonomy`, `no_match`, or `not_provided` |

`n_best_taxa` counts canonical parsed lineages; `n_weighted_best_taxa` counts full
input keys. They may differ for unsupported or differently spelled annotations
that collapse to the same conservative prefix. No biological equivalence is
inferred for mapping weights. `missing_taxonomy` takes status precedence when
both taxonomy and weight lookups are missing; both contribute to the count.

Without weights, reference/taxonomy columns duplicate the unweighted decisions,
remaining counts duplicate the original best counts, all three change/resolution
flags are false, weight measurements are blank, missing-weight count is zero
(no lookup attempted), and status is `not_provided`, including unmatched rows.
With weights but no match, remaining counts are zero, selected reference and
weight measurements are blank, taxonomy is `Unassigned`, and status is `no_match`.

## Worked synthetic example

These fixtures validate software behavior, not biological accuracy. ASV `00123`
has these matches and one ecological vector:

| Reference | Identity | Taxon | Weight |
| --- | --- | --- | --- |
| ref-a | 100% | Lactiplantibacillus plantarum | 0.05 |
| ref-b | 100% | Lactiplantibacillus plantarum | 0.05 |
| ref-c | 100% | Lactiplantibacillus paraplantarum | 0.20 |
| ref-lower | 99% | Methanobrevibacter smithii | 0.75 |

The two plantarum references share one weight of 0.05. They do not receive 0.10
combined support. The higher-weight 99% match is ineligible. Without weights,
`ref-a` is selected deterministically and taxonomy resolves only to genus. With
weights and default margins, `ref-c` wins: difference 0.15, ratio 4.

| Metadata | Unweighted comparison | Final weighted result |
| --- | --- | --- |
| Reference | ref-a | ref-c |
| Taxonomy | LCA at Lactiplantibacillus | Lactiplantibacillus paraplantarum |
| Candidate references | 3 | 1 |
| Candidate taxa | 2 | 1 |

`weight_resolved`, `assignment_changed`, and `taxonomy_changed` are true.
`min_weight_difference=0.16` or `min_weight_ratio=5` would retain all three
references and the genus LCA. Equal weights also retain the original ambiguity.
Multiple remaining references within a single species have reference ambiguity
but no disagreement between species. Deterministic reference selection remains
bookkeeping, not evidence of species/strain certainty.

The final taxonomy output for this ASV is:

```tsv
Feature ID	Taxon
00123	d__Bacteria; p__Bacillota; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactiplantibacillus; s__Lactiplantibacillus paraplantarum
```

Generated complete metadata and taxonomy examples are available in
[weighted-assignment-metadata.tsv](examples/weighted-assignment-metadata.tsv) and
[weighted-asv-taxonomy.tsv](examples/weighted-asv-taxonomy.tsv). The original
M3–M4 examples remain available to show the historical unweighted schema.

## Comparing and interpreting results

```bash
qiime tools export --input-path weighted-assignments.qza \
  --output-path weighted-assignment-export
qiime metadata tabulate --m-input-file weighted-assignments.qza \
  --o-visualization weighted-assignments.qzv
```

Compare `reference_without_weights` with `reference_with_weights`, and
`taxonomy_without_weights` with `taxonomy_with_weights`. Retain the original
unweighted ambiguity and sequence statistics when interpreting any change.
The primary taxonomy artifact always represents the final result.

Ecological priors are tie-resolution evidence, **not calibrated posterior
probabilities** or guarantees of taxonomic correctness. The method performs no
likelihood calculation, model calibration, or correction for biased training
data. Incompatible databases, incomplete annotations, or ecological priors that
poorly represent the study environment can limit interpretation. Heuristic
alignment evidence restricts weighting to best-observed matches; exhaustive
evidence is preferred. Evaluate biological accuracy separately using appropriate
held-out datasets and avoid leakage from priors into evaluation data.
