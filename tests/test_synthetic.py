"""Synthetic end-to-end tests for pmidfrag on tiny GAF + ontology files with known answers."""
import csv, gzip, json, sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from pmidfrag import __version__
from pmidfrag.cli import main

def write_ontology(path, n_terms=5):
    def node(i, lbl): return {'id': f'http://purl.obolibrary.org/obo/GO_{i:07d}', 'lbl': lbl, 'type': 'CLASS',
                              'meta': {'basicPropertyValues': [{'pred': 'http://www.geneontology.org/formats/oboInOwl#hasOBONamespace', 'val': 'biological_process'}]}}
    nodes = [node(8150, 'biological_process')] + [node(i, f'term{i}') for i in range(1, n_terms + 1)]
    edges = [{'sub': f'http://purl.obolibrary.org/obo/GO_{i:07d}', 'pred': 'is_a', 'obj': 'http://purl.obolibrary.org/obo/GO_0008150'} for i in range(1, n_terms + 1)]
    path.write_text(json.dumps({'graphs': [{'nodes': nodes, 'edges': edges}]}))

def gaf_row(prot, term, ref, ev='IDA'):
    return '\t'.join(['UniProtKB', prot, prot, 'involved_in', term, ref, ev, '', 'P', prot, '', 'protein', 'taxon:9606', '20200101', 'UniProt', '', ''])

def run_cli(tmp_path, genes, query, rows, label='out', drop_refs=(), version='2.2'):
    """Run the CLI; drop_refs removes GAF rows citing those references (independent route for checks)."""
    onto = tmp_path/'go.json'
    if not onto.exists(): write_ontology(onto)
    kept = [r for r in rows if r.split('\t')[5] not in set(drop_refs)]
    gaf = tmp_path/f'{label}.gaf'; gaf.write_text(f'!gaf-version: {version}\n' + '\n'.join(kept) + '\n')
    (tmp_path/'bg.txt').write_text('\n'.join(genes)); (tmp_path/'q.txt').write_text('\n'.join(query))
    out = tmp_path/label
    assert main(['--gaf', str(gaf), '--ontology', str(onto), '--background', str(tmp_path/'bg.txt'), '--query', str(tmp_path/'q.txt'),
                 '--out', str(out), '--min-size', '3']) == 0
    table = {r['term']: r for r in csv.DictReader(open(out/'significant_terms.csv', encoding='utf-8'))}
    with gzip.open(out/'all_terms.csv.gz', 'rt', encoding='utf-8') as f:
        all_terms = {r['term']: r for r in csv.DictReader(f)}
    return table, all_terms, out

def ids(cell):
    return set(filter(None, cell.split(';')))

def test_single_pmid_dependence(tmp_path):
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
    table, _, out = run_cli(tmp_path, genes, query, rows)
    assert 'GO:0000001' in table and 'GO:0000002' in table
    t1, t2 = table['GO:0000001'], table['GO:0000002']
    assert t1['single_pmid_dependent'] == 'True' and 'PMID:1' in ids(t1['killer_pmids'])
    assert t2['single_pmid_dependent'] == 'False'
    assert int(t1['fragility_index']) == 1 and t1['fi_status'] == 'exact'
    assert float(t1['evidence_concentration']) == pytest.approx(5/6)
    # supporting PMIDs are written as a list, and the count matches it
    assert ids(t1['support_pmids']) == {'PMID:1', 'PMID:2'} and int(t1['n_support_pmids']) == 2
    assert ids(t2['support_pmids']) == {f'PMID:{100+j}' for j in range(6)} and int(t2['n_support_pmids']) == 6
    assert int(t1['hits_without_pmid']) == 0
    # term2 needs at least two deletions; the status must match the value
    if t2['fragility_index']:
        fi = int(t2['fragility_index']); assert fi >= 2
        assert t2['fi_status'] == ('exact' if fi == 2 else 'upper_bound')
    else:
        assert t2['fi_status'].startswith('not_abolished')
    for t in table.values():  # killer split is a partition; direct killers always support the term
        assert ids(t['killer_pmids']) == ids(t['killer_pmids_direct']) | ids(t['killer_pmids_via_bh'])
        assert not ids(t['killer_pmids_direct']) & ids(t['killer_pmids_via_bh'])
        assert ids(t['killer_pmids_direct']) <= ids(t['support_pmids'])
    summ = json.loads((out/'summary.json').read_text())
    assert summ['significant_terms'] == 2 and summ['single_pmid_dependent_terms'] == 1
    assert summ['version'] == __version__ and 'interpretation' not in summ
    assert (out/'report.html').exists() and (out/'ledger.jsonl.gz').exists()
    assert 'replicat' not in (out/'report.html').read_text(encoding='utf-8')

def test_dependence_through_bh_only(tmp_path):
    """X (term2) and Z (term3) keep each other significant through BH: deleting a PMID that supports only Z
    abolishes X although X's own counts do not change, and vice versa."""
    genes = [f'G{i}' for i in range(60)]; query = genes[:10]; rows = []
    for i in range(7): rows.append(gaf_row(genes[i], 'GO:0000001', f'PMID:{10 + i}'))   # A: K=8, k=7
    rows.append(gaf_row('G20', 'GO:0000001', 'PMID:17'))
    for j, g in enumerate(['G7', 'G8', 'G9']): rows.append(gaf_row(g, 'GO:0000002', f'PMID:{30 + j}'))  # X: K=5, k=3
    for g in ['G21', 'G22']: rows.append(gaf_row(g, 'GO:0000002', 'PMID:33'))
    for j, g in enumerate(['G0', 'G1', 'G2', 'G3']): rows.append(gaf_row(g, 'GO:0000003', f'PMID:{40 + j}'))  # Z: K=9, k=4
    for g in genes[23:28]: rows.append(gaf_row(g, 'GO:0000003', 'PMID:44'))
    for g in genes[30:40]: rows.append(gaf_row(g, 'GO:0000004', 'PMID:50'))  # filler, not significant
    table, base_all, out = run_cli(tmp_path, genes, query, rows)
    X, Z = table['GO:0000002'], table['GO:0000003']
    x_own, z_own = {'PMID:30', 'PMID:31', 'PMID:32'}, {'PMID:40', 'PMID:41', 'PMID:42', 'PMID:43'}
    assert ids(X['killer_pmids_direct']) == x_own and ids(X['killer_pmids_via_bh']) == z_own
    assert ids(Z['killer_pmids_direct']) == z_own and ids(Z['killer_pmids_via_bh']) == x_own
    assert ids(table['GO:0000001']['killer_pmids']) == set()
    summ = json.loads((out/'summary.json').read_text())
    assert summ['single_pmid_dependent_terms'] == 2 and summ['single_pmid_dependent_terms_via_bh_only'] == 0
    # independent route: drop the GAF rows of PMID:40 and rerun; X keeps K and k but is no longer significant
    _, after, _ = run_cli(tmp_path, genes, query, rows, label='drop40', drop_refs=['PMID:40'])
    xb, xa = base_all['GO:0000002'], after['GO:0000002']
    assert (xa['K'], xa['k']) == (xb['K'], xb['k']) and float(xa['q']) > 0.05 >= float(xb['q'])

def test_fi_status_and_hits_without_pmid(tmp_path):
    genes = [f'G{i}' for i in range(80)]; query = genes[:15]; rows = []
    for i in range(6):  # R: every query hit supported by two PMIDs, one IEA-only hit, IEA-only background
        rows.append(gaf_row(genes[i], 'GO:0000001', f'PMID:{60 + i}')); rows.append(gaf_row(genes[i], 'GO:0000001', f'PMID:{70 + i}'))
    rows.append(gaf_row('G6', 'GO:0000001', 'GO_REF:0000002', 'IEA'))
    for g in ['G40', 'G41']: rows.append(gaf_row(g, 'GO:0000001', 'GO_REF:0000002', 'IEA'))
    for i in range(12): rows.append(gaf_row(genes[i], 'GO:0000002', f'PMID:{100 + i}'))  # S: survives five deletions
    for g in ['G42', 'G43']: rows.append(gaf_row(g, 'GO:0000002', 'PMID:120'))
    for g in ['G0', 'G1', 'G2']: rows.append(gaf_row(g, 'GO:0000003', 'PMID:200'))  # T: needs exactly two deletions
    for g in ['G3', 'G4', 'G5']: rows.append(gaf_row(g, 'GO:0000003', 'PMID:201'))
    rows.append(gaf_row('G6', 'GO:0000003', 'PMID:202'))
    for g in genes[44:48]: rows.append(gaf_row(g, 'GO:0000003', 'PMID:203'))
    for g in genes[50:60]: rows.append(gaf_row(g, 'GO:0000004', 'PMID:300'))  # filler
    table, _, out = run_cli(tmp_path, genes, query, rows)
    R, S, T = table['GO:0000001'], table['GO:0000002'], table['GO:0000003']
    assert R['fi_status'] == 'not_abolished_shortlist_exhausted' and R['fragility_index'] == '' and R['greedy_removal_set'] == ''
    assert int(R['hits_without_pmid']) == 1 and int(R['n_support_pmids']) == 12 and R['single_pmid_dependent'] == 'False'
    assert S['fi_status'] == 'not_abolished_after_5' and S['fragility_index'] == '' and len(ids(S['greedy_removal_set'])) == 5
    assert T['fi_status'] == 'exact' and int(T['fragility_index']) == 2 and ids(T['greedy_removal_set']) == {'PMID:200', 'PMID:201'}
    summ = json.loads((out/'summary.json').read_text())
    assert summ['fi_status_counts'] == {'exact': 1, 'not_abolished_after_5': 1, 'not_abolished_shortlist_exhausted': 1}
    # independent route: removing both greedy PMIDs from the GAF abolishes T; removing either alone does not
    for drop, lost in ((['PMID:200', 'PMID:201'], True), (['PMID:200'], False), (['PMID:201'], False)):
        _, after, _ = run_cli(tmp_path, genes, query, rows, label='drop_' + '_'.join(d[5:] for d in drop), drop_refs=drop)
        assert (float(after['GO:0000003']['q']) > 0.05) == lost

def without_qualifier(rows):
    return ['\t'.join(z[:3] + [''] + z[4:]) for z in (r.split('\t') for r in rows)]

def test_gaf21_empty_qualifier_and_input_errors(tmp_path):
    """GAF 2.1 leaves the qualifier of positive annotations empty: the same rows must give the same result as
    GAF 2.2 rows with involved_in. Negative control: an empty qualifier in a GAF 2.2 file is not accepted.
    A query with no background gene and an empty term family exit with code 2 and a message, not a traceback."""
    genes = [f'G{i}' for i in range(60)]; query = genes[:10]; rows = []
    for g in query[:5]: rows.append(gaf_row(g, 'GO:0000001', 'PMID:1'))
    rows.append(gaf_row(query[5], 'GO:0000001', 'PMID:2'))
    for g in genes[20:24]: rows.append(gaf_row(g, 'GO:0000001', 'PMID:3'))
    for j, g in enumerate(query[:6]): rows.append(gaf_row(g, 'GO:0000002', f'PMID:{100+j}'))
    for g in genes[30:34]: rows.append(gaf_row(g, 'GO:0000002', 'PMID:200'))
    for g in genes[50:56]: rows.append(gaf_row(g, 'GO:0000004', 'PMID:301'))
    t22, a22, o22 = run_cli(tmp_path, genes, query, rows, label='v22')
    t21, a21, o21 = run_cli(tmp_path, genes, query, without_qualifier(rows), label='v21', version='2.1')
    assert t21.keys() == t22.keys() and len(t22) >= 1
    for term in t22:
        for col in ('q', 'k', 'K', 'support_pmids', 'killer_pmids', 'fragility_index', 'fi_status'):
            assert t21[term][col] == t22[term][col], (term, col)
    assert a21 == a22
    qc21 = json.loads((o21/'summary.json').read_text())['ledger_qc']
    qc22 = json.loads((o22/'summary.json').read_text())['ledger_qc']
    assert qc21['gaf_version'] == '2.1' and qc21['empty_qualifier_read_as'] == 'involved_in'
    assert qc22['gaf_version'] == '2.2' and qc22['empty_qualifier_read_as'] is None
    assert qc21['accepted_rows'] == qc22['accepted_rows'] == len(rows)
    onto = tmp_path/'go.json'
    (tmp_path/'bg.txt').write_text('\n'.join(genes)); (tmp_path/'q.txt').write_text('\n'.join(query))
    (tmp_path/'none.txt').write_text('X1\nX2\n')
    bad = tmp_path/'bad22.gaf'; bad.write_text('!gaf-version: 2.2\n' + '\n'.join(without_qualifier(rows)) + '\n')
    def rc(gaf, q, label):
        return main(['--gaf', str(gaf), '--ontology', str(onto), '--background', str(tmp_path/'bg.txt'),
                     '--query', str(q), '--out', str(tmp_path/label), '--min-size', '3'])
    assert rc(bad, tmp_path/'q.txt', 'bad22') == 2                       # empty family: nothing accepted
    assert rc(tmp_path/'v22.gaf', tmp_path/'none.txt', 'noquery') == 2   # no query gene in the background
    summ = json.loads((o22/'summary.json').read_text())['inputs']
    assert summ['propagation'] == 'is_a+part_of' and summ['entity'] == 'protein'
    assert len(summ['background_sha256']) == len(summ['query_sha256']) == 64 and summ['mapping_sha256'] is None
