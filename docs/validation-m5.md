# M5 implementation and validation

Recorded on 2026-11-08. M5 extends the existing assignment action with optional
single-vector ecological priors. No enhanced OTU aggregation (M6), complete
reassignment workflow (M7), or combined pipeline (M8) was implemented.

## Implementation

- `_weighting.py` uses the standard library to validate full-lineage keys,
  reject normalization collisions, compare one weight per exact-best taxon,
  apply both margins, preserve missing evidence, and materialize final decisions.
- `_assignment.py` retains the unweighted resolver and adds an optional weighted
  pass over exact-best hits. Integer-ratio identity is rechecked. The returned
  mapping and taxonomy represent the final decision.
- `_assign.py` reads a single standard BIOM vector and serializes the existing
  SequenceCharacteristics and Taxonomy formats. No new type or dependency.
- `plugin_setup.py` adds optional `FeatureTable[RelativeFrequency]` input and
  three parameters to `assign-closed-reference`. No new action is registered.
- `tests/test_weighting.py` adds 22 tests using existing synthetic evidence and
  small reusable weight fixtures.
- The assignment/weighting guides, README, changelog, and generated example
  TSVs describe the complete signature, rules, and 43-column schema.

The original 30 metadata columns preserve their values. Thirteen columns are
added in all outputs. Without priors, assignments and remaining counts are
explicitly equivalent, weight metrics are missing, and status is `not_provided`.
One M3 test asserting that weighted columns do not yet exist was updated to
check these specified M5 semantics. All other existing tests remain intact.
The conventional clustering and alignment algorithms are unchanged.

## Executed results

| Check | Environment | Result |
| --- | --- | --- |
| Complete repository regression suite, including the original upstream URL-based usage example | QIIME 2 2024.10, Python 3.10.14, VSEARCH 2.22.1 | **143 passed**, 140 warnings |
| Enhanced alignment, assignment, and weighting suites | Rachis/QIIME 2 2026.7, Python 3.12.13 | **73 passed**, 79 subtests passed, 3 warnings |
| Assignment and weighting suites from a fresh wheel, outside the repository | QIIME 2 2024.10 | **51 passed**, 126 warnings |
| Original unweighted results compared with installed pre-M5 implementation | Standard-library core, 64 configurations | Exact equality of all original summary fields, taxonomy, and mapping |
| Fresh-wheel artifact smoke test and generated documentation examples | QIIME 2 2024.10 | Passed; 43 columns, all six original ASVs, unchanged evidence fingerprints, external commands prohibited during assignment |
| Isolated CLI cache refresh and assignment help | QIIME 2 2024.10 | Passed; optional RelativeFrequency input and margin/policy parameters exposed |
| Weighted assignment CLI, VSEARCH absent from PATH | QIIME 2 2024.10 | Passed; both output artifacts saved |
| Standard-library core without QIIME 2 in executing interpreter | System Python 3.13 | Passed |
| Wheel build/install, flake8, compilation, whitespace checks | Development workspace | Passed |

Existing conda environments were not modified;
the wheel installation and CLI cache were isolated under temporary directories.

No failures remain in the final automated selections. The first sandboxed test
attempt could not import QIIME 2 because system CPU queries were denied; the
same selection was rerun outside that sandbox and passed. Existing
QIIME/pandas deprecation warnings, intentionally empty FASTA warnings, and
single-rank Unassigned taxonomy warnings remain visible. No warning suppression
was used to obtain these results.

The 64 pre-M5 comparisons cover eight taxonomy inputs, both selection policies,
and four identity thresholds (0.90, 0.97, 0.995, 1.0), including missing and
incomplete taxonomy and numeric-looking identifiers. The comparison projects
out the 13 newly requested columns and checks every original field.

## Coverage and invariants

Synthetic tests cover unique references, multiple references in the winning
taxon, taxonomic ambiguity, higher-weight lower-identity/near-best references,
uniform priors, all-zero candidates, zero runners-up, tiny/subnormal positive
weights and ratio overflow, both inclusive margin boundaries, no per-ASV
renormalization, missing weights/taxonomy, incompatible database lineages,
whitespace aliases and collisions, reference multiplicity/order, incomplete
lineages, deterministic/excluded references, threshold invariance from 90–100%,
unmatched ASVs, ID preservation, taxonomy/metadata artifact validation and
round trips, and input vector shape validation.

The exact-best set is recomputed using integer ratios before a weight can be
looked up. Reference multiplicity never adds support. Missing candidates veto
ranking rather than disappearing or receiving zeros. Ties and failed margins
retain all exact-best candidate taxa. The unweighted classifications and all
original ambiguity fields remain available alongside final weighted taxonomy.

Every assignment/weighting test prohibits external command execution using a
subprocess guard after framework initialization. Evidence files are hashed
before and after weighted/unweighted calls. Assignment neither realigns nor
modifies the stored evidence. The schema-1 evidence and original fixtures were
not rewritten.

## Reproducible commands

In a dedicated QIIME 2 environment containing the original q2-vsearch:

```bash
python -m pip install . --no-deps
qiime dev refresh-cache
qiime meta-vsearch assign-closed-reference --help
python -m pytest q2_meta_vsearch/tests/test_weighting.py -x --tb=short -q
python -m pytest q2_meta_vsearch/tests --tb=short -q
```

The full suite uses the preserved upstream example's existing QIIME tutorial
URLs. A local-only selection can omit that example:

```bash
python -m pytest q2_meta_vsearch/tests \
  --ignore=q2_meta_vsearch/tests/test_examples.py --tb=short -q
```

The [weighting guide](weighting.md) contains the exact action signature, CLI
example, complete additional metadata schema, and explanation of missing values.
Its generated example files come from actual synthetic action outputs, without
presenting those fixtures as biological performance evidence.

## Limitations and next milestone

Patched QIIME 2 2026.10 and QIIME 2 2025.10 remain untested. The 2026.7 environment
was tested for enhanced actions, not the entire legacy plugin. A real researcher
provided Clawback artifact and biological accuracy dataset have not been tested.
Compatibility uses the standard RelativeFrequency type and a synthetic vector;
no q2-clawback import or runtime requirement was introduced.

Taxonomy mapping normalizes whitespace only. Different prefixes, database
versions, synonyms, missing ranks, and trailing delimiters are not automatically
reconciled. Conservative taxonomy parsing may report an ancestor despite a
single weighted input key. Unknown taxonomies cannot be assigned invented ranks.
Weighting memory is proportional to the exact-best references/taxa for one query,
plus the taxonomy mapping, weight vector, and output summaries. Other hits remain
in disk-backed SQLite. Priors are tie-resolution evidence, not calibrated
posterior probabilities.

After explicit authorization, M6 should consume the final internal ASV-to-reference
mapping and reuse conventional BIOM aggregation. Validate per-sample read-count
conservation, deterministic/excluded assignments, representative sequences, and
reference-level taxonomy while retaining the ASV audit metadata. M6 is not part
of this change.
