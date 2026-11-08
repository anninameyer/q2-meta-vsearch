# q2-meta-vsearch

A standalone QIIME 2 plugin forked from
[q2-vsearch](https://github.com/qiime2/q2-vsearch). The distribution is
`q2-meta-vsearch`, its Python package is `q2_meta_vsearch`, and its CLI namespace
is `qiime meta-vsearch`. It can be installed alongside the original plugin.

The enhanced workflow is being developed incrementally as **ALIGN → ASSIGN →
CLUSTER**. Ecological priors will resolve exact sequence ties only; they must
never override stronger sequence evidence. Ecological weighting, enhanced OTU
clustering, ecological reassignment, and the combined workflow are future work.

M1–M4 are implemented: standalone plugin, reusable exhaustive alignment,
ASV assignment metadata, and unweighted taxonomy-aware assignment.

```bash
qiime meta-vsearch align-closed-reference \
  --i-sequences asvs.qza \
  --i-reference-sequences references.qza \
  --p-identity-floor 0.90 \
  --p-search-mode exhaustive \
  --p-strand plus \
  --p-threads 1 \
  --o-alignments alignments.qza \
  --o-blast6 alignments-blast6.qza
```

All accepted matches are retained, including lower-identity matches. Exact-best
flags use integer alignment counts, never rounded BLAST6 identity. Every ASV
appears in the query inventory, including unmatched ASVs. Heuristic mode is
available for speed and does not guarantee globally optimal matches. See the
[evidence schema and search contract](docs/alignment-evidence.md).

Assign references and LCA taxonomy from that evidence without running VSEARCH:

```bash
qiime meta-vsearch assign-closed-reference \
  --i-alignments alignments.qza \
  --i-reference-taxonomy reference-taxonomy.qza \
  --p-identity-threshold 0.97 \
  --p-identity-tolerance 0 \
  --o-assignment-metadata assignments.qza \
  --o-taxonomy asv-taxonomy.qza
```

Exact-best references alone determine taxonomy. Default reference selection is
lexicographic bookkeeping; ambiguity remains explicit. See the
[assignment guide, complete metadata schema, and worked example](docs/assignment.md)
and [M3–M4 validation report](docs/validation-m3-m4.md).

The conventional actions, including `cluster-features-closed-reference`, retain
the original q2-vsearch algorithms. Their representative sequences are input
ASVs relabeled with reference IDs. Chimera statistics use the independently
namespaced `MetaVsearchUchimeStats` type and formats to avoid collisions with
q2-vsearch.

## Installation and validation

Install into a dedicated QIIME 2 environment, keeping the distribution-provided
q2-vsearch installed:

```bash
python -m pip install . --no-deps
qiime dev refresh-cache
qiime meta-vsearch --help
qiime vsearch --help
python -m pytest q2_meta_vsearch/tests --tb=short
```

Dependencies are managed by `conda-recipe/meta.yaml` and the QIIME 2 distribution.
The primary target is patched QIIME 2 2026.10; older-release compatibility is
claimed only when explicitly validated. See `docs/validation.md` for reproducible
validation instructions and recorded limitations.
