# pmidfrag

Provenance audits for Gene Ontology over-representation analysis. For every enriched term, `pmidfrag` reports
which publications the result stands on:

| Output column | Meaning |
|---|---|
| `support_pmids` | the distinct PMIDs that support the query hits of the term after ontology propagation (`;`-separated); `n_support_pmids` is their number |
| `hits_without_pmid` | number of query hits of the term supported only by annotations that cite no PMID (e.g. `GO_REF` rows); no publication deletion can remove them |
| `single_pmid_dependent` | removing all annotations citing **one** PMID (anywhere in the background) makes the term non-significant |
| `killer_pmids` | the PMIDs with that property, split into `killer_pmids_direct` (the deletion removes at least one query hit of the term) and `killer_pmids_via_bh` (the deletion removes none of the term's query hits; the term is lost through the BH re-adjustment of the family) |
| `fragility_index`, `fi_status` | number of PMIDs whose joint removal makes the term non-significant, from a greedy search (below) |
| `greedy_removal_set` | the PMIDs removed by the greedy search |
| `evidence_concentration` | largest share of the term's query hits removed by deleting a single PMID |
| `hits_only_IEP_HEP` | query hits supported only by expression-pattern evidence (IEP/HEP) |
| `q_without_IEP_HEP` | BH-adjusted p after dropping IEP/HEP annotations |

All deletions act on the whole background, closures are recomputed and the full fixed term family is re-adjusted with BH,
so the numbers are exact for the declared annotation table, not approximations.

**Fragility index.** FI 1 is exact (single-publication dependence is tested for every candidate PMID). For other
terms a greedy search removes, one at a time, the PMID that raises the term's q most, choosing among the (at most 30)
PMIDs whose deletion on its own changes the term's counts. `fi_status` is `exact` for FI 1 and 2, `upper_bound` for
FI 3-5, `not_abolished_after_5` when five removals leave the term significant, and
`not_abolished_shortlist_exhausted` when fewer than five such PMIDs exist. The last two are not bounds: PMIDs that
change the term's counts only jointly (hits supported by two or more publications) are not searched.

**What it is not.** Dependence on one publication is a pointer to what to re-read; it is not an error flag. In the
evaluations of the accompanying article, dependent terms were no more likely than terms with similar q values to fail
replication in independent cohorts or to be lost after annotation updates. Use it to decide where to look.

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
- `--mapping`: `accession<TAB>gene` TSV (e.g. UniProt -> Ensembl/Entrez). Accessions that the table maps to more
  than one gene, or to a gene outside the background, are dropped. Omit when the gene ids are the GAF object ids
  themselves (e.g. `--db SGD` with SGD ids).
- `--db`, `--taxon`, `--entity`: GAF filters (defaults `UniProtKB`, `9606`, `protein`).
- Projection: biological process, positive qualifiers (`involved_in` and `acts_upstream_of*`), `is_a` + `part_of`
  propagation, NOT/ND/extension rows excluded, family = terms with 5-500 background hits. Change with `--relations`,
  `--no-part-of`, `--min-size`, `--max-size`.
- GAF versions: 2.2, and 2.0/2.1, whose positive annotations have an empty qualifier; there an empty qualifier is
  read as `involved_in` (the version is taken from the `!gaf-version` header and recorded in `summary.json`).
- The run stops with exit code 2 and a message if no query gene is in the background or no term is in the family.

Outputs: `significant_terms.csv`, `all_terms.csv.gz`, `summary.json` (with input SHA-256 and QC ledger), `report.html`,
`ledger.jsonl.gz` (every GAF row with its inclusion or exclusion reason, evidence code and references).

## Reference run

`examples/airway/` contains the airway (GSE52778) dexamethasone gene lists, the identifier mapping and the audited
per-term results; `bash examples/airway/run_airway.sh` downloads the pinned GO 2026-06-19 files, verifies their
SHA-256 values, runs `pmidfrag` and checks the output against the reference table (see `examples/airway/README.md`).

## Tests

```bash
python -m pytest -q tests/              # synthetic end-to-end tests with known answers
bash examples/airway/run_airway.sh      # real-data regression (downloads ~83 MB)
```

## Changes

- 0.1.2: GAF 2.0/2.1 input is read (0.1.1 rejected their empty qualifiers and stopped with an error); exit code 2
  with a message for a query outside the background or an empty term family; `summary.json` also records the
  SHA-256 values of the background, query and mapping files, the entity and propagation settings and the GAF
  version. Results for GAF 2.2 input are unchanged (airway regression).
- 0.1.1: per-term supporting-PMID lists and hits without PMID support; killer PMIDs split into direct and BH-only;
  `fi_status` for every fragility index; FI threshold follows `--alpha`; more tests; study-specific helpers removed.
  Significance, dependence flags, killer PMIDs and fragility-index values are unchanged from 0.1.0 at the default
  `--alpha`; the former `>5` label is now an empty value with its `fi_status`.
- 0.1.0: first release.

## Citation

If you use pmidfrag, please cite the software (see `CITATION.cff`) and the accompanying article (in preparation).

## License

MIT for the code (see `LICENSE`). GO annotations and ontology are CC BY 4.0 (Gene Ontology Consortium); cite the GO release you use.
