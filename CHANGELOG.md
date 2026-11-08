# q2-meta-vsearch changelog

## Unreleased — M3 and M4

- Add unweighted `assign-closed-reference` using existing evidence and reference
  taxonomy, without VSEARCH execution or a frequency-table input.
- Report one metadata-compatible summary per original ASV, exact identity
  counts, reference/taxonomic ambiguity, LCA, diagnostics, and deterministic
  reference selection or exclusion.
- Add conservative taxonomy parsing and separate species-alternative/ssg.
  labels. Near-best tolerance is diagnostic only.
- Preserve numeric-looking identifiers through canonical artifact metadata views
  using additive directory-format transformers for existing semantic types.
- Add synthetic evidence, assignment invariants, metadata round-trip tests,
  worked examples, and documentation. No ecological weighting or enhanced OTU
  aggregation has been added.

## Unreleased — M1 and M2

- Namespace the standalone distribution, Python package, plugin, and chimera
  statistics formats/types so q2-vsearch can remain installed alongside it.
- Preserve conventional clustering, dereplication, merging, chimera filtering,
  and sequence-statistics algorithms.
- Add `align-closed-reference`, defaulting to exhaustive VSEARCH search and
  retaining all accepted hits above a configurable 90% default identity floor.
- Add versioned reusable evidence with integer statistics, exact-best flags,
  complete query inventories, input fingerprints, and search settings.
- Add standard BLAST6 projection, synthetic integration tests, evidence
  validation, installation instructions, and a synthetic profiling helper.
- Taxonomy, ecological weighting, reassignment, and enhanced clustering remain
  future milestones.
