# Reference dataset: airway, dexamethasone up-regulated genes

Inputs derived from the Bioconductor `airway` package 1.32.0 (GSE52778; four donor cell lines, dexamethasone
vs. control), DESeq2 1.52.0 with design `~ cell + dex`, BH <= 0.05 and |log2FC| >= 1:

| File | Content |
|---|---|
| `background.txt` | 12,595 measured genes (Ensembl) with an unambiguous UniProt mapping |
| `query_up.txt`, `query_down.txt` | 433 up- and 410 down-regulated genes |
| `uniprot_to_ensembl.tsv` | `UniProt<TAB>Ensembl` mapping used for the background |
| `expected_terms_up.csv.gz` | audited per-term results for the up list (all 6,462 terms of the family) |

Run (Linux/macOS shell with `curl` and `sha256sum`, after `pip install .` in the repository root):

```bash
bash examples/airway/run_airway.sh
```

The script downloads `HUMAN-uniprot.gaf.gz` and `go-basic.json` from the GO 2026-06-19 release, verifies their
SHA-256 values, runs `pmidfrag` and compares the output with `expected_terms_up.csv.gz`. Expected summary:
family 6,462 terms; 1,070 candidate PMIDs; 18 significant terms; 8 single-publication dependent (1 only through the
BH re-adjustment); 2 of the top 10 dependent; no PMID abolishes all terms; no term lost without IEP/HEP; FI status
10 exact, 1 upper bound, 1 not abolished after five removals, 6 with the shortlist exhausted. The comparison checks K, k, q and q without
IEP/HEP for every term (tolerance 1e-12) and dependence flags, killer PMIDs, support counts and ranks for the
significant terms, plus the internal consistency of the supporting-PMID lists, the killer split and the FI status. `python examples/airway/check_regression.py <out> --tamper` must report FAIL (negative control).
