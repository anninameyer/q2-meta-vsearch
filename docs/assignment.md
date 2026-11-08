# Unweighted assignment: M3 and M4

`assign-closed-reference` consumes existing alignment evidence and reference
taxonomy. It returns ASV metadata and hierarchical ASV taxonomy. It accepts no
feature table or ecological weights and never invokes VSEARCH. Enhanced OTU
clustering and ecological reassignment workflows remain out of scope.

## Action signature

```text
assign_closed_reference(
    alignments: FeatureData[ReferenceAlignments],
    reference_taxonomy: FeatureData[Taxonomy],
    identity_threshold: Float = 0.90,
    identity_tolerance: Float = 0.0,
    unresolved_policy: Choices['deterministic-reference', 'exclude']
        = 'deterministic-reference',
    max_explicit_species: Int = 2
) -> (
    assignment_metadata: FeatureData[SequenceCharacteristics],
    taxonomy: FeatureData[Taxonomy]
)
```

`identity_threshold` is a fraction in `(0, 1]` and must be at least the stored
alignment identity floor. `identity_tolerance` is in **percentage points**, from
0 to 100. `max_explicit_species` must be a positive integer.

```bash
qiime meta-vsearch assign-closed-reference \
  --i-alignments alignments.qza \
  --i-reference-taxonomy reference-taxonomy.qza \
  --p-identity-threshold 0.97 \
  --p-identity-tolerance 0.2 \
  --p-unresolved-policy deterministic-reference \
  --p-max-explicit-species 2 \
  --o-assignment-metadata assignments.qza \
  --o-taxonomy asv-taxonomy.qza
```

## Review of alignment evidence

Schema-1 M1–M2 evidence supplies everything required for these milestones:

1. Exact-best flags use integer match/alignment counts. Assignment independently
   recomputes the maximum as an exact rational; rounded identity is never used
   to decide ties, thresholds, or near-best eligibility.
2. `hits.tsv` retains every genuine accepted match above the search floor, not
   just top hits. Assignment does not modify this file.
3. `queries.tsv` inventories every original ASV, including queries without hits.
4. Standard BLAST6 no-hit records are outside `hits.tsv`. The new action accepts
   only `FeatureData[ReferenceAlignments]`, and rejects BLAST6 artifacts. Generic
   external BLAST6 consumers must recognize no-hit sentinels (target `*`, zero
   alignment length); no universal behavior is claimed for external software.
5. For fixed evidence, raising a threshold either retains the same maximum and
   exact-best set, or makes the query unmatched. It cannot promote a lower hit.
   Exhaustive evidence provides this guarantee relative to the eligible
   reference database; heuristic evidence provides it only for observed hits.

No M1–M2 alignment algorithm or evidence schema change was necessary. The new
action validates the complete artifact before assignment. The standard-library
core checks individual records, inventory fingerprints, and hit membership/counts
when used directly, and recomputes decisions rather than trusting stored flags.

## Decision rules

Among alignments meeting the assignment threshold, find the maximum exact
`matches/alignment_length`. Only references with that exact identity enter
reference selection or taxonomic consensus. Each reference counts once even
when both strands match; alignment counts remain separate.

The selected reference is the lexicographically first original reference ID.
Reference FASTA ordering, query abundance, and reference multiplicity do not
break taxonomic ties. Selection is an OTU bookkeeping decision, not evidence
that the ASV uniquely identifies the selected species or strain.

Taxonomic consensus is the **lowest common ancestor (LCA)**: the longest shared,
contiguous, explicitly annotated lineage prefix across distinct exact-best
candidate lineages. It is not majority voting. Ten references annotated as A
and one as B give the same LCA as one A and one B.

Reference ambiguity and taxonomic ambiguity differ: multiple exact-best
references with the same complete species lineage have reference ambiguity,
but support an unambiguous species annotation.

With `exclude`, references are not selected if known candidate lineages differ
or any exact-best reference lacks a usable lineage. Counts, taxonomy, and
sequence statistics remain present. Identical incomplete annotations support
their shared stated prefix; incompleteness alone is not disagreement, and does
not trigger exclusion. It never implies resolution at an unannotated species.

## Near-best diagnostics

For a best identity of 0.998 (99.8%) and tolerance 0.2 percentage points,
diagnostic candidates extend down to 0.996 (99.6%), subject to the assignment
threshold and stored search floor. `n_near_best_reference_hits` counts only
**additional, strictly lower-identity references**, excluding every reference
that already has an exact-best alignment. With zero tolerance it is zero.

Changing tolerance cannot change reference selection, LCA, or future
weight-eligible candidates. A threshold of 0.997 would clip these diagnostics
at 99.7%. Hits not retained by the original alignment run cannot be recovered
by setting a larger tolerance.

## Taxonomy validation and conservative interpretation

Taxonomy identifiers are read as strings. Duplicate/empty IDs and malformed rows
are rejected. Additional taxonomy IDs are allowed: reference databases normally
contain many records without accepted hits. Missing exact-best reference IDs,
empty annotations, and annotations without a usable root are counted explicitly.

If any exact-best candidate has no usable lineage, the LCA is `Unassigned`:
the unknown candidate could belong to any lineage. The chosen reference remains
available under the default policy, with a missing-taxonomy status. Missing
taxonomy never causes a lower-identity annotated reference to replace a best hit.
Missing annotations on lower-identity hits do not influence consensus.

Supported conventions are `d__/p__/c__/o__/f__/g__/s__`, explicit kingdom roots
such as `k__Fungi`, legacy `k__Bacteria` and other recognized domain roots,
`D_0__` through `D_6__`, and unprefixed
seven-rank paths beginning with a recognized domain. Output is normalized to
`d__/p__/c__/o__/f__/g__/s__` (or `k__` for an explicit kingdom root); labels
remain case-sensitive. A common `k__Fungi` ancestor is reported as Kingdom, never
inferred as Domain. Prefixes identify ranks, not spelling guesses. An isolated
binomial, out-of-order ranks, and unsupported conventions do not get guessed
ranks. Use explicitly formatted
reference taxonomy when those conventions are ambiguous.

Consensus stops at a missing, placeholder, or unsupported rank. A gap in the
middle is not filled using names at later ranks. Thus an explicit genus after
an omitted phylum does not manufacture a complete lineage. This conservative
policy can report a higher ancestor than a less strict interpretation.

Species reporting requires an explicit genus and a conservative binomial or
species epithet at the species rank. A bare epithet is expanded using that
explicit genus. `sp.`, unknown/uncultured annotations, strain/organism labels,
numeric GTDB species placeholders, and unsupported name forms are not promoted
to named species. Complex names such as some Candidatus annotations may therefore
remain at genus. No synonym or database-version reconciliation is attempted.

### Human-readable labels

- One supported species: `Lactiplantibacillus plantarum`.
- Two species in the same genus:
  `Lactiplantibacillus paraplantarum|plantarum`.
- More species than `max_explicit_species` in the same genus:
  `Lactiplantibacillus ssg.`.
- Across genera, or with unknown species: the name of the lowest shared rank.

Species epithets are sorted lexicographically for stable labels. `ssg.` means
unresolved among multiple named species in one genus; it is **not a formally
defined taxonomic subgroup**. These labels are diagnostic metadata only.
Hierarchical `FeatureData[Taxonomy]` always contains the LCA, never `ssg.` or a
species-alternative expression.

## Complete metadata schema

Every original ASV appears once. Categorical booleans use `true`/`false`, as
QIIME 2 metadata supports categorical and numeric columns rather than boolean
columns. Unavailable numbers and references are blank; counts are zero for
unmatched ASVs. Taxonomy uses the explicit `Unassigned` label.

| Column | Type / units | Meaning |
| --- | --- | --- |
| `matched` | categorical boolean | At least one stored hit meets the assignment threshold |
| `best_identity` | numeric fraction | Exact maximum, serialized as a convenient floating value |
| `best_identity_numerator` | integer-valued numeric | Raw matching columns for the statistics representative |
| `best_identity_denominator` | integer-valued numeric | Raw internal alignment length for that representative; ratio need not be reduced |
| `n_accepted_reference_hits` | integer-valued numeric | Distinct references meeting the assignment threshold |
| `n_best_reference_hits` | integer-valued numeric | Distinct exact-best reference IDs |
| `n_best_alignments` | integer-valued numeric | Exact-best alignment rows, counting strands separately |
| `n_near_best_reference_hits` | integer-valued numeric | Additional lower-identity reference IDs in the diagnostic band |
| `best_query_coverage` | numeric fraction | Paired nucleotide columns / query length |
| `best_reference_coverage` | numeric fraction | Paired nucleotide columns / reference length |
| `best_alignment_length` | integer-valued numeric | Alignment columns excluding terminal gaps |
| `best_mismatches` | integer-valued numeric | Representative mismatch count |
| `best_gaps` | integer-valued numeric | Representative gap-base count excluding terminal gaps |
| `best_strand` | categorical | Representative query orientation, `+` or `-` |
| `reference_without_weights` | categorical | Selected reference, or blank when unmatched/excluded |
| `n_remaining_reference_hits` | integer-valued numeric | Exact-best reference candidates; no pruning occurs in M3–M4, even under exclusion |
| `reference_assignment_status` | categorical | `unique_reference`, `deterministic_reference`, `excluded_taxonomic_ambiguity`, or `no_match` |
| `n_best_taxa` | integer-valued numeric | Distinct known canonical candidate lineage prefixes; missing annotations are not invented taxa |
| `taxonomy_without_weights` | categorical | Unweighted hierarchical LCA |
| `taxonomy_lca` | categorical | Same LCA, explicitly named for auditing |
| `taxonomy_ambiguity_label` | categorical | Human-readable species alternatives or ancestor label |
| `taxonomic_assignment_status` | categorical | `resolved`, `incomplete`, `ambiguous`, `missing_taxonomy`, `partial_missing_taxonomy`, or `no_match` |
| `taxonomic_resolution_rank` | categorical | Species, Genus, Family, Order, Class, Phylum, Kingdom (explicit root only), Domain, or Unassigned |
| `search_mode` | categorical | Original exhaustive or heuristic search |
| `search_complete` | categorical boolean | Original search-completeness flag; heuristic results remain best-observed |
| `identity_threshold` | numeric fraction | Applied assignment threshold |
| `identity_tolerance` | numeric percentage points | Diagnostic tolerance |
| `reference_selection_policy` | categorical | Requested `deterministic-reference` or `exclude` policy |
| `statistics_reference` | categorical | Reference supplying the representative alignment statistics |
| `n_missing_taxonomy_reference_hits` | integer-valued numeric | Exact-best reference IDs with absent/empty/unusable lineage annotations |

For every `best_*` alignment statistic, use the lexicographically first exact-best
reference and, if both orientations tie for it, the plus strand first. This does
not favor a reference with higher coverage or a longer alignment. Under exclusion,
the same diagnostic representative remains named in `statistics_reference`
even though `reference_without_weights` is blank. The original evidence contains
all alternative statistics. Unmatched ASVs have no representative statistics.

## Worked synthetic example

These labels and sequences are software-validation fixtures, not biological
accuracy results. ASV `00123` has three exact-best references:

| Reference | Identity | Species |
| --- | --- | --- |
| ref-a | 100% | Lactiplantibacillus plantarum |
| ref-b | 100% | Lactiplantibacillus plantarum |
| ref-c | 100% | Lactiplantibacillus paraplantarum |

Its metadata reports `n_best_reference_hits=3`, `n_best_alignments=3`,
`n_best_taxa=2`, selected reference `ref-a`, taxonomic status `ambiguous`, and
resolution `Genus`. The ambiguity label is
`Lactiplantibacillus paraplantarum|plantarum`. Reference multiplicity does not
give plantarum two votes. The taxonomy output is:

```tsv
Feature ID	Taxon
00123	d__Bacteria; p__Bacillota; c__Bacilli; o__Lactobacillales; f__Lactobacillaceae; g__Lactiplantibacillus
```

The fixture also has an accepted 99% hit annotated as Archaea. It contributes to
the accepted-hit count and, with sufficient tolerance, the near-best count, but
never affects the selected reference or this LCA.

Actual generated full examples are in [examples/assignment-metadata.tsv](examples/assignment-metadata.tsv)
and [examples/asv-taxonomy.tsv](examples/asv-taxonomy.tsv).

## Export and metadata operations

```bash
qiime tools export --input-path assignments.qza --output-path assignment-export
qiime tools export --input-path asv-taxonomy.qza --output-path taxonomy-export
qiime metadata tabulate --m-input-file assignments.qza \
  --o-visualization assignments.qzv
```

The summary file is `assignment-export/sequence_characteristics.tsv`. It is
standard typed QIIME 2 metadata, including a `#q2:types` directive (not an ASV
row), so numeric-looking selected reference IDs stay categorical. Read it with
`qiime2.Metadata.load`, or use the artifact's metadata view. Generic CSV readers
must recognize metadata directives and preserve ID columns as strings.

The plugin adds direct directory-format views to Metadata and DataFrame for the
existing SequenceCharacteristics and canonical taxonomy formats, plus a safe
taxonomy Series view. These preserve
IDs such as `00123` and `NA`, and honor metadata type directives. They do not
replace q2-types' file-level transformers or register competing formats/types.
Existing file-level transformer paths outside the artifact view remain owned by
q2-types; prefer the artifact metadata view or `Metadata.load` for this summary.

## Internal architecture and later work

`_assignment.py` and `_taxonomy.py` use ordinary Python objects and the standard
library. The core returns an `AssignmentResult` containing summary records,
ASV taxonomy, and an ASV-to-reference dictionary (selected ASVs only). `_assign.py`
provides the QIIME 2 adapter and serialization. The mapping is represented in
the summary's reference column rather than a new mapping semantic type.

Evidence is indexed in temporary disk-backed SQLite; query hits are streamed
twice and not accumulated across ASVs. The taxonomy mapping and final ASV summary
are held in memory. Assignment validates evidence but never realigns it.

The same evidence can be passed to this action with another supported threshold,
tolerance, policy, or taxonomy without modifying the evidence. A lower threshold
than the stored floor is rejected. This is evidence reuse within the M3–M4 action;
no ecological reassignment workflow or feature-table aggregation is implemented.

M5 should insert an optional taxon-level decision after exact-best candidate
identification and before final consensus/reference materialization, reusing
these same threshold, diagnostic, and taxonomy functions. A lower-identity hit
must never become eligible. The schema has no placeholder weight columns; the
requested weighted columns can be added when M5 is authorized.
