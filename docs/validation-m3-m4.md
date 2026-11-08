# M3–M4 implementation and validation

Recorded on 2026-11-08. Implemented ASV assignment metadata and unweighted
taxonomy-aware assignment only. No M5 weights, M6 enhanced OTU aggregation, or
M7 ecological reassignment workflow has been implemented. Conventional clustering
remains available unchanged.

## Evidence review and preserved implementation

The review confirmed that schema-1 evidence supplies unrounded integer counts,
all accepted hits above the recorded floor, and every query ID. Assignment reads
only genuine `hits.tsv` rows, rejects BLAST6 artifacts, and independently
recomputes exact-best candidates. A higher threshold either preserves an ASV's
exact-best set or makes it unmatched; it cannot replace that set with lower hits.

No alignment algorithm or evidence schema correction was required. The M1–M2
alignment modules, evidence validators/formats/types, legacy algorithm modules,
and existing alignment/clustering regression tests were compared byte for byte
against the previously installed M1–M2 wheel and remain unchanged.

The important output-view issue identified during review was pandas inference
of numeric-looking IDs and NA spellings. Additive canonical directory-format
transformers now preserve IDs in metadata, DataFrame, and taxonomy Series views.
They do not replace existing file-level transformer registrations. Typed metadata
directives also preserve numeric-looking selected reference IDs as categorical.

## Implemented files and interface

- `_assignment.py`: standard-library evidence reader and exact-candidate resolver;
  returns summary records, ASV taxonomy, and a selected ASV-to-reference mapping.
- `_taxonomy.py`: conservative lineage parsing, distinct-lineage LCA, resolution
  ranks, and human-readable ambiguity labels, with no voting.
- `_assign.py`: QIIME 2 adapter and metadata/taxonomy serialization.
- `plugin_setup.py`: registers `assign-closed-reference` using existing semantic
  types and no frequency-table input.
- `_transformer.py`: identifier-preserving views of existing canonical formats.
- `tests/test_assign.py` and `tests/data/assignment/`: synthetic evidence,
  taxonomies, ordinary-object core tests, action/serialization checks, and
  strict no-external-command guards.
- `docs/assignment.md`, `docs/examples/`, README, changelog, and the historical
  evidence/validation links: complete schema, worked examples, usage, and limits.

The exact action signature and all 30 metadata columns are documented in
[assignment.md](assignment.md). The two extra audit columns beyond the requested
28 are `statistics_reference` and `n_missing_taxonomy_reference_hits`. No weighted
placeholder columns were added. The internal mapping is available through the
core result and the summary reference column; it does not need a new semantic type.

## Executed results

| Check | Environment | Result |
| --- | --- | --- |
| Entire repository suite, including the preserved upstream URL-based usage example | QIIME 2 2024.10, Python 3.10.14, VSEARCH 2.22.1 | **121 passed** |
| Enhanced alignment and assignment suites | Rachis/QIIME 2 2026.7, Python 3.12.13, VSEARCH 2.22.1 | **51 passed**, 49 subtests passed |
| Assignment suite from an installed wheel, outside the repository | QIIME 2 2024.10 | **29 passed** |
| Real assignment CLI with VSEARCH absent from PATH | QIIME 2 2024.10 | Passed |
| Standard CLI export and `qiime metadata tabulate` on the output artifact, with VSEARCH absent from PATH | QIIME 2 2024.10 | Passed |
| Standard-library core without QIIME 2 installed in the executing interpreter | System Python 3.13 | Passed |
| Wheel contents, flake8, compilation, whitespace checks | Development workspace | Passed |

No failures remain in these selections. During development, a test needed a
string path for the older Metadata.save API, and the newer framework needed
platform-provenance initialization before the external-command guard. These
were test-only compatibility corrections. The guard still prohibits every
external command during the assignment checks.

Observed warnings come from existing QIIME/pandas APIs, intentionally empty
alignment fixtures, and the newer taxonomy validator warning about single-rank
`Unassigned` annotations. They are not hidden or converted into passing claims.

Existing environments were not modified. Wheel installations, CLI caches, and
integration artifacts used temporary directories. The synthetic evidence fixtures
were generated using real VSEARCH before assignment tests; assignment tests do
not run VSEARCH. The existing upstream usage test was also run successfully and
is no longer excluded from this milestone's final full-suite result.

Patched QIIME 2 2026.10 and QIIME 2 2025.10 remain unavailable locally and untested.
The 2026.7 minimal environment was tested for enhanced actions, not all legacy
actions. VSEARCH 2.30.4 has not been executed. No compatibility or biological
performance claim is made for untested environments or real microbial datasets.

## Coverage and invariants

The assignment tests cover unique matches; same-species reference ambiguity;
two/three species; shared genus/family/domain; cross-domain ambiguity; explicit
fungal kingdoms; missing species/ranks/annotations; absent reference IDs;
nonstandard species labels; reference multiplicity; FASTA permutations;
deterministic/excluded references; near-best diagnostics; exact precision hidden
by BLAST6; thresholds 90–100% and an exact 99.5% boundary; unmatched queries;
numeric/NA identifiers; representative coverage; metadata export/reimport and
artifact round trips; taxonomy validation; and BLAST6 semantic rejection.

Every assignment test installs a subprocess guard after framework initialization.
Tests assert that no external command is invoked. The CLI test additionally
removed VSEARCH from PATH while executing assignment, export, and metadata
tabulation. The core contains no QIIME 2, pandas, or subprocess imports.

## Reproducible commands

In a dedicated QIIME 2 environment containing the original q2-vsearch:

```bash
python -m pip install . --no-deps
qiime dev refresh-cache
qiime meta-vsearch assign-closed-reference --help

python -m pytest q2_meta_vsearch/tests/test_assign.py -x --tb=short -q
python -m pytest q2_meta_vsearch/tests --tb=short -q
```

The final command includes the upstream usage example and therefore needs network
access to its existing QIIME tutorial URLs. For a strictly local selection:

```bash
python -m pytest q2_meta_vsearch/tests \
  --ignore=q2_meta_vsearch/tests/test_examples.py --tb=short -q
```

Target-environment creation instructions remain in
[the M1–M2 report](validation.md#target-installation-and-tests). The provided
2026.10 Linux environment file is still untested locally.

## Methodological limits and next steps

Taxonomy is normalized conservatively. Unknown ranks are not inferred from
species names; unsupported annotations can reduce resolution. Identical partial
lineages support only their explicit common prefix. Missing best-hit annotations
make the LCA Unassigned. The reference taxonomy's database version cannot be
verified from labels alone; this action checks identifiers, not the biological
correctness of an annotation.

Near-best counts include additional references only, clipped to the assignment
threshold and original evidence floor. `n_remaining_reference_hits` remains the
exact-best candidate count even when exclusion suppresses selection. Explicit
kingdom roots can produce a Kingdom rank instead of falsely inferring Domain.
These conventions and representative-statistic selection are documented in the
schema and action help.

The next authorized milestone would be M5: accept one standard RelativeFrequency
taxon-weight vector, validate/normalize full taxonomy keys and missing weights,
compare each exact-best taxon once, and add the weighted audit columns while
preserving unweighted outputs. Extend tests with uniform and extreme priors,
weight margins/ratios, and the invariant that no lower-identity or diagnostic
near-best hit can be selected. Do not introduce OTU aggregation in that milestone.
