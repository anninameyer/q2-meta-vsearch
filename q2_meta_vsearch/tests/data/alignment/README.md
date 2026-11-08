# Synthetic alignment fixtures

These are software-validation data, not biological observations or evidence of
real-world taxonomic performance. Sequences were generated with Python's
`random.Random(812)` (200-base core followed by 3,000-base precision sequence).
The unrelated sequence uses `random.Random(97513)`.

- `queries` / `references`: two identical references plus a lower-identity
  reference; includes a query with no accepted match.
- `suboptimal`: two dispersed substitutions in the better reference versus
  three clustered substitutions in the worse reference. Shared-word ranking
  selects the worse reference first under conventional stopping at 90% identity.
  The individual-reference files serve as a per-reference search oracle.
- `precision`: two versus three substitutions across 3,000 bases. Both BLAST6
  identities print as 99.9%, but their exact ratios differ.
- `geometry`: forward/reverse complements, a 160-base subsequence, a 14-base
  internal insertion, and a query containing N.
- `reference-ambiguous`: all Ns; demonstrates k-mer-filter bypass and VSEARCH's
  ambiguous-character identity behavior, not a useful biological reference.
- `identifiers`: identical sequences with identifiers that must not be coerced
  to numbers, missing values, or VSEARCH abundance annotations.
- `duplicate-ids` / `empty`: invalid-input fixtures.
