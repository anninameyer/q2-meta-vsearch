# q2-meta-vsearch changelog

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
