# Alignment evidence, schema 1

`align-closed-reference` is the enhanced functionality currently implemented.
It accepts ASV sequences and references, without a feature table, taxonomy, or
weights. Assignment and enhanced clustering are not implemented in M1–M2.

```text
align_closed_reference(
    sequences: FeatureData[Sequence],
    reference_sequences: FeatureData[Sequence],
    identity_floor: Float = 0.90,
    strand: Choices['plus', 'both'] = 'plus',
    search_mode: Choices['exhaustive', 'heuristic'] = 'exhaustive',
    threads: Threads = 1
) -> (
    alignments: FeatureData[ReferenceAlignments],
    blast6: FeatureData[BLAST6]
)
```

Identity floor is a fraction in `(0, 1]`; thread count 0 uses available cores.

## Search behavior

The action reuses the original `run_command` infrastructure and uses VSEARCH
`--usearch_global`. Exhaustive mode explicitly sets:

```text
--maxaccepts 0 --maxrejects 0 --minwordmatches 0 --maxhits 0
```

No `--top_hits_only` filter is used. Every accepted alignment at or above the
identity floor is retained. Exhaustive means every input reference is eligible
for pairwise alignment on every configured query strand, under VSEARCH's
scoring model; it does not enumerate alternative score-optimal alignments for
the same pair. VSEARCH performs one alignment per reference and strand.

Heuristic mode uses `maxaccepts=1`, `maxrejects=32`, and `minwordmatches=12`.
It retains the hits that this search accepts, but **cannot guarantee retrieval
of globally best references or all equally best references**. The original
conventional clustering action is also preserved unchanged.

Both modes explicitly preserve the conventional alignment model:

| Setting | Value |
| --- | --- |
| Identity definition | 2 |
| Match / mismatch | +2 / -4 |
| Gap opening / extension | `20I/2E` / `2I/1E` |
| Query / database masking | none / DUST |
| Word length | 8 |
| Minimum query/reference coverage | 0 / 0 |
| Minimum sequence length | 1 |

Maximum sequence length is set to the longest input sequence, so VSEARCH's
default 50,000-base limit cannot silently drop longer input records. Empty
sequences/collections and duplicate identifiers within a collection are rejected.

The installed VSEARCH log supplies the version. Versions before 2.3.1 are
rejected because the zero-word bypass was not reliable there; this is a
capability floor, **not a claim that every later version has been tested**.
The installed 2.22.1 manual and synthetic integration tests were used to verify
the zero stopping limits, reporting limits, and k-mer bypass. See
[VSEARCH documentation](https://torognes.github.io/vsearch/commands/vsearch-usearch_global.1.html).

## Evidence directory

`FeatureData[ReferenceAlignments]` uses a namespaced
`MetaVsearchReferenceAlignmentsDirFmt` containing exactly three files.

### `hits.tsv`

One row per accepted query/reference/strand alignment, sorted by original query
ID, reference ID, and strand. No-hit queries produce no rows here. Columns:

```text
query_id reference_id strand matches alignment_length mismatches
gap_bases gap_opens gap_extensions query_length reference_length
query_start query_end reference_start reference_end raw_score cigar
identity query_coverage reference_coverage is_exact_best
```

The actual delimiter is a tab and the file has a header.

- `alignment_length` excludes terminal gaps, following VSEARCH's userout field.
- `matches`, `mismatches`, and gap counts are integers.
- `identity = matches / alignment_length`, stored as a fraction, not a percent.
- Coverage fractions are `(matches + mismatches) / sequence_length`.
- `gap_bases`, openings, and extensions exclude terminal gaps.
- Coordinates are VSEARCH's trimmed `qilo/qihi/tilo/tihi` fields, relative to the
  oriented sequences. The query orientation is explicit in `strand`.
- `cigar` is the original VSEARCH `caln` field, including terminal gaps. Its
  I/D convention is VSEARCH's convention, not a substitute for SAM CIGAR.
- `raw_score` is the VSEARCH alignment score, not a BLAST bit score.

`is_exact_best` is 1 exactly when the integer match/length ratio equals the
maximum for that query among retained hits. Comparisons use integer cross
multiplication. With heuristic evidence, this means **best observed**, and does
not assert a global optimum. Both strands can contribute exact-best rows.

Printed identity and coverage columns are conveniences. The integer counts are
authoritative. The pure `identity_fraction` helper returns an exact rational.

### `queries.tsv`

Exactly one row per original ASV, sorted by original ID:

```text
query_id sequence_length sequence_sha256 accepted_alignments status
```

`status` is `matched` or `no_hit`. There is no taxonomy status at this stage.
`accepted_alignments` counts alignment rows, not distinct reference IDs.

FASTA IDs are interpreted using QIIME 2/scikit-bio's ID convention; FASTA
descriptions are not part of a feature ID. Both queries and references are
relabelled with temporary, reversible safe tokens before VSEARCH. Identifiers
such as `00123`, `NA`, and `asv;size=7;odd` are restored without parsing them as
numbers or VSEARCH annotations. Duplicate sequences with different IDs remain
separate records.

### `search.json`

Records the evidence schema version, installed VSEARCH version, search mode,
completeness, identity definition/floor, strand, effective options/flags,
retention policy, sanitized command/log, and elapsed VSEARCH wall time.

Each input collection has its record count, maximum length, order-independent
SHA256 fingerprint, and an input-order fingerprint. The order fingerprint is
particularly relevant to heuristic candidate ordering.

The fingerprint algorithm is versioned: uppercase each sequence and compute its
ASCII SHA256; sort records by original identifier; serialize each
`[identifier, sequence_length, sequence_sha256]` as compact ASCII-escaped JSON
followed by LF; SHA256 the concatenated bytes. The order fingerprint uses input
order instead of sorted ID order. Thus sequence wrapping and case do not change
the fingerprint, but identifier or nucleotide changes do.

Full format validation checks schema/settings consistency, fingerprints, query
membership/lengths, hit counts, unique query/reference/strand rows, integer
statistics, CIGAR consistency, derived values, and exact-best flags. Parsing,
sorting, and cross-file validation use temporary disk-backed SQLite storage.
Python does not accumulate all reference sequences or alignment rows in memory.
VSEARCH itself can still require substantial memory and computation.

## Standard BLAST6 projection

The second output uses the existing q2-types `BLAST6Format`. Its twelve columns
are projected from the exact evidence using VSEARCH's BLAST6 conventions:
one-decimal identity percentages, internal alignment length, mismatch/gap-opening
counts, global coordinates, E-value -1, and bit score 0.

Unmatched queries have explicit no-hit records with target `*` and zero alignment
statistics, as with VSEARCH `--output_no_hits`. They are not fabricated
alignments and never appear in `hits.tsv`. This also permits a valid BLAST6
artifact when every query is unmatched. The evidence artifact is authoritative.

## Constraints for later milestones

- Priors may resolve **exact-best** sequence ties only.
- Current identity tolerance is zero. A future nonzero tolerance may label
  near-best matches for reporting, without making them weight-eligible.
- A complete 90% run can be filtered to higher identity thresholds without
  rerunning VSEARCH. Reuse below the recorded floor is invalid.
- Changing references, strands, scoring, or identity definition requires new
  alignments. Changing taxonomy or ecological priors does not.
- Reference-level and taxonomic ambiguity are separate. Future deterministic
  reference selection defaults to stable reference-ID ordering and retains
  unresolved taxonomy; it is OTU bookkeeping, not strain identification.
- Low coverage, terminal-gap exclusion, and ambiguous-base matching can yield
  high identity with weak biological specificity. These are preserved VSEARCH
  behaviors, explicitly exposed in the evidence, not biological validation.
