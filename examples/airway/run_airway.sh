#!/usr/bin/env bash
# Reference run: airway (GSE52778) dexamethasone up-regulated genes against the GO 2026-06-19 release.
# Usage (from the repository root, after `pip install .`): bash examples/airway/run_airway.sh
# Downloads two GO files (~83 MB), verifies their SHA-256, runs pmidfrag and compares with the audited table.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
GO_DIR="${GO_DIR:-$HERE/go_2026-06-19}"
OUT="${OUT:-$HERE/out_up}"
BASE=https://release.geneontology.org/2026-06-19
mkdir -p "$GO_DIR"

fetch() {  # url file sha256
  if [ ! -s "$GO_DIR/$2" ]; then
    echo "[airway] downloading $2"
    curl -sSfL --retry 3 -o "$GO_DIR/$2" "$1"
  fi
  echo "$3  $GO_DIR/$2" | sha256sum -c -
}
fetch "$BASE/annotations/gaf/HUMAN-uniprot.gaf.gz" HUMAN-uniprot.gaf.gz 258f6ea163375c036929478c400f25bf32e6d83df547cc4b864dfb4ec20c1ffb
fetch "$BASE/ontology/go-basic.json" go-basic.json a336667a9582a5230b54d854ab83163c64488455f39bb2e7cefbbdca821b9d3a

echo "[airway] running pmidfrag (about 20-40 s)"
pmidfrag --gaf "$GO_DIR/HUMAN-uniprot.gaf.gz" --ontology "$GO_DIR/go-basic.json" \
  --background "$HERE/background.txt" --query "$HERE/query_up.txt" \
  --mapping "$HERE/uniprot_to_ensembl.tsv" --out "$OUT"

echo "[airway] comparing with the audited reference table"
python "$HERE/check_regression.py" "$OUT"
