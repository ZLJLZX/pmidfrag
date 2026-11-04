# pmidfrag

Provenance audits for Gene Ontology over-representation analysis. For every enriched term, `pmidfrag` reports
which publications the result stands on:

| Output column | Meaning |
|---|---|
| `support_pmids` | distinct PMIDs that support the query hits of the term (after ontology propagation) |
| `single_pmid_dependent` | removing all annotations citing **one** PMID (anywhere in the background) makes the term non-significant |
| `killer_pmids` | the PMIDs with that property |
| `fragility_index` | minimum number of PMIDs whose removal makes the term non-significant (greedy, exact for 1, `>5` if more than 5) |
| `evidence_concentration` | largest share of the term's query hits attributable to a single PMID |
| `hits_only_IEP_HEP` | query hits supported only by expression-pattern evidence (IEP/HEP) |
| `q_without_IEP_HEP` | BH-adjusted p after dropping IEP/HEP annotations |

All deletions act on the whole background, closures are recomputed and the full fixed term family is re-adjusted with BH,
so the numbers are exact for the declared annotation table, not approximations.

**What it is not.** Dependence on one publication is a pointer to what to re-read; it is not an error flag. In our
evaluations it does not predict whether a term replicates in an independent dataset (replication was governed by the
q-value margin), and expression-evidence circularity was absent in the test case. Use it to decide where to look.

## Install

```bash
pip install .            # Python >= 3.11, numpy, scipy
# or
docker build -t pmidfrag .
```

## Run

```bash
pmidfrag --gaf HUMAN-uniprot.gaf.gz --ontology go-basic.json \
         --background measured_genes.txt --query up_genes.txt \
         --mapping uniprot_to_gene.tsv --out results/
```

- `--background`: one gene id per line, the measured genes (never the whole genome).
- `--query`: gene ids, subset of the background.
- `--mapping`: `accession<TAB>gene` TSV (e.g. UniProt -> Ensembl/Entrez). Accessions mapping to more than one
  background gene are dropped. Omit when the gene ids are the GAF object ids themselves (e.g. `--db SGD` with SGD ids).
- `--db`, `--taxon`, `--entity`: GAF filters (defaults `UniProtKB`, `9606`, `protein`).
- Projection: biological process, positive qualifiers (`involved_in` and `acts_upstream_of*`), `is_a` + `part_of`
  propagation, NOT/ND/extension rows excluded, family = terms with 5-500 background hits. Change with `--relations`,
  `--no-part-of`, `--min-size`, `--max-size`.

Outputs: `significant_terms.csv`, `all_terms.csv.gz`, `summary.json` (with input SHA-256 and QC ledger), `report.html`,
`ledger.jsonl.gz` (every accepted GAF row with evidence code and references).

## Reference run

`examples/airway/` contains the airway (GSE52778) dexamethasone gene lists, the identifier mapping and the audited
per-term results; `bash examples/airway/run_airway.sh` downloads the pinned GO 2026-06-19 files, verifies their
SHA-256 values, runs `pmidfrag` and checks the output against the reference table (see `examples/airway/README.md`).

## Tests

```bash
python -m pytest -q tests/              # synthetic end-to-end test
bash examples/airway/run_airway.sh      # real-data regression (downloads ~83 MB)
```

## Citation

If you use pmidfrag, please cite the software (see `CITATION.cff`) and the accompanying article (in preparation).

## License

MIT for the code (see `LICENSE`). GO annotations and ontology are CC BY 4.0 (Gene Ontology Consortium); cite the GO release you use.
