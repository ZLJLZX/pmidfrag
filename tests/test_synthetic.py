"""Synthetic end-to-end test for pmidfrag: a tiny GAF + ontology where one term depends on one PMID."""
import gzip, json, csv, sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from pmidfrag.cli import main

def write_ontology(path):
    def node(i, lbl): return {'id': f'http://purl.obolibrary.org/obo/GO_{i:07d}', 'lbl': lbl, 'type': 'CLASS',
                              'meta': {'basicPropertyValues': [{'pred': 'http://www.geneontology.org/formats/oboInOwl#hasOBONamespace', 'val': 'biological_process'}]}}
    nodes = [node(8150, 'biological_process')] + [node(i, f'term{i}') for i in range(1, 6)]
    edges = [{'sub': f'http://purl.obolibrary.org/obo/GO_{i:07d}', 'pred': 'is_a', 'obj': 'http://purl.obolibrary.org/obo/GO_0008150'} for i in range(1, 6)]
    path.write_text(json.dumps({'graphs': [{'nodes': nodes, 'edges': edges}]}))

def gaf_row(prot, term, ref, ev='IDA'):
    return '\t'.join(['UniProtKB', prot, prot, 'involved_in', term, ref, ev, '', 'P', prot, '', 'protein', 'taxon:9606', '20200101', 'UniProt', '', ''])

def test_single_pmid_dependence(tmp_path):
    onto = tmp_path/'go.json'; write_ontology(onto)
    N = 60; genes = [f'G{i}' for i in range(N)]; query = genes[:10]
    rows = []
    # term1: 6 of 10 query genes annotated, 5 via PMID:1 only, 1 via PMID:2 ; 4 background genes outside query
    for g in query[:5]: rows.append(gaf_row(g, 'GO:0000001', 'PMID:1'))
    rows.append(gaf_row(query[5], 'GO:0000001', 'PMID:2'))
    for g in genes[20:24]: rows.append(gaf_row(g, 'GO:0000001', 'PMID:3'))
    # term2: 6 query genes each with its own PMID -> robust
    for j, g in enumerate(query[:6]): rows.append(gaf_row(g, 'GO:0000002', f'PMID:{100+j}'))
    for g in genes[30:34]: rows.append(gaf_row(g, 'GO:0000002', 'PMID:200'))
    # filler terms so the family has >= 3 terms with 5..500 hits
    for g in genes[40:50]: rows.append(gaf_row(g, 'GO:0000003', 'PMID:300', 'IEP'))
    for g in genes[50:56]: rows.append(gaf_row(g, 'GO:0000004', 'PMID:301'))
    gaf = tmp_path/'test.gaf'; gaf.write_text('!gaf-version: 2.2\n' + '\n'.join(rows) + '\n')
    (tmp_path/'bg.txt').write_text('\n'.join(genes)); (tmp_path/'q.txt').write_text('\n'.join(query))
    out = tmp_path/'out'
    assert main(['--gaf', str(gaf), '--ontology', str(onto), '--background', str(tmp_path/'bg.txt'), '--query', str(tmp_path/'q.txt'),
                 '--out', str(out), '--min-size', '3']) == 0
    table = {r['term']: r for r in csv.DictReader(open(out/'significant_terms.csv'))}
    assert 'GO:0000001' in table and 'GO:0000002' in table
    assert table['GO:0000001']['single_pmid_dependent'] == 'True' and 'PMID:1' in table['GO:0000001']['killer_pmids']
    assert table['GO:0000002']['single_pmid_dependent'] == 'False'
    assert int(table['GO:0000001']['fragility_index']) == 1 and int(table['GO:0000002']['fragility_index']) >= 2
    assert float(table['GO:0000001']['evidence_concentration']) == pytest.approx(5/6)
    summ = json.loads((out/'summary.json').read_text()); assert summ['significant_terms'] == 2 and summ['single_pmid_dependent_terms'] == 1
    assert (out/'report.html').exists() and (out/'ledger.jsonl.gz').exists()
