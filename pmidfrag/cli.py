"""pmidfrag: single-publication dependence and evidence-type audits for GO over-representation analysis.

Inputs: a GAF (2.0-2.2; an empty 2.0/2.1 qualifier is read as involved_in), a GO JSON ontology (go-basic.json[.gz]), a background gene list, a query gene list,
and (unless genes are the GAF object identifiers themselves) a two-column mapping TSV `accession<TAB>gene`.
Outputs: per-term CSV with q, the supporting PMIDs, query hits without PMID support, single-publication
dependence (killer PMIDs split into direct and BH-only), fragility index with its status, evidence
concentration and q after dropping IEP/HEP; JSON summary; a self-contained HTML report.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, html, json, re, sys, time
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from . import __version__
from .core import load_ontology, closure, hit_counts, bh, sha256
from .ledger import build_ledger, EXPRESSION_CODES
from .engine import FragilityEngine, pmids

ISO = re.compile(r'-\d+$')
FI_STEPS = 5
FIELDS = ['rank', 'term', 'name', 'q', 'p', 'K', 'k', 'n_support_pmids', 'support_pmids', 'hits_without_pmid',
          'single_pmid_dependent', 'killer_pmids', 'killer_pmids_direct', 'killer_pmids_via_bh',
          'fragility_index', 'fi_status', 'greedy_removal_set', 'evidence_concentration',
          'hits_only_IEP_HEP', 'q_without_IEP_HEP']
NOTE = ('Each row is a significant term. Deleting a PMID removes every annotation that cites it from the whole '
        'background, recomputes propagation and re-applies BH to the fixed term family. Single-publication '
        'dependence marks the publications a term rests on; it does not by itself show that an annotation is wrong.')

def read_list(p):
    return [l.strip() for l in Path(p).read_text(encoding='utf-8').splitlines() if l.strip() and not l.startswith('#')]

def load_p2g(mapping, background):
    """accession -> gene, keeping only accessions that the table maps to exactly one gene, if that gene is in the background."""
    bg = set(background)
    if mapping is None:
        return {g: g for g in bg}
    genes_of = defaultdict(set)
    with open(mapping, encoding='utf-8') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            acc, gene = line.rstrip('\n').split('\t')[:2]
            genes_of[ISO.sub('', acc.strip().removeprefix('UniProtKB:'))].add(gene.strip())
    return {a: next(iter(gs)) for a, gs in genes_of.items() if len(gs) == 1 and next(iter(gs)) in bg}

def fragility_index(engine, query, base_q, i, shortlist, max_steps=FI_STEPS, alpha=0.05):
    """Greedy search for the fewest PMIDs whose joint deletion makes term i non-significant.

    Called only for terms that no single deletion abolishes. Each step deletes the shortlisted PMID that raises
    q_i most. Returns (fi, removed, status): 'exact' for fi == 2 (single deletions are known not to suffice),
    'upper_bound' for fi 3..max_steps, and fi None with 'not_abolished_after_<max_steps>' or
    'not_abolished_shortlist_exhausted' (fewer than max_steps shortlisted PMIDs). Neither of the last two is a
    bound: PMIDs that change the term's counts only jointly (redundant support) are not shortlisted.
    """
    removed = []; cache = {}
    def q_after(rem):
        key = tuple(sorted(rem))
        if key in cache: return cache[key]
        rs = set(rem); q = engine.query_set(query); qset = set(q)
        affected = set().union(*(engine.sources.get(s, set()) for s in rs)) if rs else set()
        base = engine.baseline(q); K = base['K'].copy(); k = base['k'].copy()
        for g in affected:
            kept = engine._close(r for r in engine.rows[g] if not (rs & set(r['references'])))
            gone = engine.closed[g] - kept
            if gone:
                idx = list(gone); K[idx] -= 1
                if g in qset: k[idx] -= 1
        cache[key] = bh(engine._tail(len(q), K, k)); return cache[key]
    for _ in range(max_steps):
        best = None
        for s in shortlist:
            if s in removed: continue
            r = q_after(removed + [s])
            if best is None or r[i] > best[1][i]: best = (s, r)
        if best is None:
            return None, removed, 'not_abolished_shortlist_exhausted'
        removed.append(best[0])
        if best[1][i] > alpha:
            return len(removed), removed, ('exact' if len(removed) <= 2 else 'upper_bound')
    return None, removed, 'not_abolished_after_%d' % max_steps

def run(args):
    t0 = time.perf_counter(); out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    ont = load_ontology(args.ontology, use_part_of=not args.no_part_of)
    background = read_list(args.background); query_in = read_list(args.query)
    p2g = load_p2g(args.mapping, background)
    bg = {g for g in p2g.values()}
    ledger_path = out/'ledger.jsonl.gz'
    if ledger_path.exists(): ledger_path.unlink()
    rows, ledger_qc = build_ledger(args.gaf, p2g, ont, args.relations.split(','), ledger_path,
                                   db=args.db, taxon=args.taxon, entity=args.entity)
    byg = defaultdict(list)
    for r in rows: byg[r['gene']].append(r['term'])
    closed = {g: closure(byg.get(g, []), ont.ancestors) for g in bg}
    allterms = sorted({t for s in closed.values() for t in s}); ti = {t: i for i, t in enumerate(allterms)}
    K = hit_counts(closed, bg, ti)
    family = [t for t in allterms if t != 'GO:0008150' and args.min_size <= K[ti[t]] <= args.max_size]
    query = sorted(set(query_in) & bg); unmapped = sorted(set(query_in) - bg)
    if not family:
        print('pmidfrag: no GO term has %d-%d background hits (%d accepted GAF rows; row reasons: %s)'
              % (args.min_size, args.max_size, len(rows), json.dumps(ledger_qc['reason_counts'])), file=sys.stderr)
        return 2
    if not query:
        print('pmidfrag: none of the %d query genes is in the background (%d background genes after mapping)'
              % (len(query_in), len(bg)), file=sys.stderr)
        return 2
    engine = FragilityEngine(rows, bg, ont.ancestors, family)
    res = engine.evaluate(query, args.alpha, source_sink=lambda *_: None)
    base = res['baseline']; sig = res['baseline_significant_indices']; qs = engine.query_set(query)
    sup = res['support']
    # evidence concentration and expression-evidence exclusion
    maxlost = np.zeros(len(family), dtype=np.int64)
    for s in res['candidate_pmids']:
        lost, _, _ = engine.source_delta(s); cnt = np.zeros(len(family), dtype=np.int64)
        for g in qs:
            if g in lost: cnt[lost[g]] += 1
        np.maximum(maxlost, cnt, out=maxlost)
    expr_ids = [r['row_id'] for r in rows if r['evidence'] in EXPRESSION_CODES]
    from .engine import delta_for_row_ids, score_delta
    E = score_delta(engine, query, base, delta_for_row_ids(engine, expr_ids))
    table = []
    for rank, i in enumerate(sig, 1):
        per = []
        for s in res['candidate_pmids']:
            lost, ix, _ = engine.source_delta(s)
            if i in set(ix.tolist()):
                per.append((sum(1 for g in qs if g in lost and i in set(lost[g].tolist())), s))
        per.sort(reverse=True); shortlist = [s for _, s in per[:30]]
        if res['fragile1'][i]:
            fi, removed, status = 1, res['killers'][i][:1], 'exact'
        else:
            fi, removed, status = fragility_index(engine, query, base['q'], i, shortlist, alpha=args.alpha)
        table.append({'rank': rank, 'term': family[i], 'name': ont.names.get(family[i], family[i]),
                      'q': float(base['q'][i]), 'p': float(base['p'][i]), 'K': int(base['K'][i]), 'k': int(base['k'][i]),
                      'n_support_pmids': int(sup['support_pmids'][i]), 'support_pmids': ';'.join(sup['support_pmid_ids'][i]),
                      'hits_without_pmid': int(sup['only_no_pmid_hit_genes'][i]),
                      'single_pmid_dependent': bool(res['fragile1'][i]), 'killer_pmids': ';'.join(res['killers'][i]),
                      'killer_pmids_direct': ';'.join(res['killers_direct'][i]),
                      'killer_pmids_via_bh': ';'.join(res['killers_via_bh'][i]),
                      'fragility_index': '' if fi is None else fi, 'fi_status': status,
                      'greedy_removal_set': ';'.join(removed),
                      'evidence_concentration': float(maxlost[i]/base['k'][i]) if base['k'][i] else 0.0,
                      'hits_only_IEP_HEP': int(sup['only_IEP_HEP_hit_genes'][i]),
                      'q_without_IEP_HEP': float(E['q'][i])})
    with open(out/'significant_terms.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); w.writerows(table)
    with gzip.open(out/'all_terms.csv.gz', 'wt', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['term', 'name', 'K', 'k', 'p', 'q', 'q_without_IEP_HEP'])
        for i, t in enumerate(family): w.writerow([t, ont.names.get(t, t), int(base['K'][i]), int(base['k'][i]), float(base['p'][i]), float(base['q'][i]), float(E['q'][i])])
    summary = {'tool': 'pmidfrag', 'version': __version__,
               'inputs': {'gaf': str(args.gaf), 'gaf_sha256': sha256(args.gaf), 'ontology': str(args.ontology),
               'ontology_sha256': sha256(args.ontology), 'background_sha256': sha256(args.background),
               'query_sha256': sha256(args.query), 'mapping': args.mapping,
               'mapping_sha256': sha256(args.mapping) if args.mapping else None,
               'background_genes': len(background), 'background_used': len(bg),
               'query_genes': len(query_in), 'query_used': len(query), 'query_unmapped': unmapped,
               'relations': args.relations.split(','), 'db': args.db, 'taxon': args.taxon, 'entity': args.entity,
               'propagation': 'is_a' if args.no_part_of else 'is_a+part_of',
               'family_size_bounds': [args.min_size, args.max_size], 'alpha': args.alpha},
               'ledger_qc': ledger_qc, 'family_size': len(family), 'candidate_pmids': res['candidate_count'],
               'query_rows_without_pmid': res['query_accepted_rows_without_pmid'],
               'significant_terms': len(sig), 'single_pmid_dependent_terms': int(res['fragile1_count']),
               'single_pmid_dependent_terms_via_bh_only': sum(1 for i in sig if res['killers'][i] and not res['killers_direct'][i]),
               'top10_single_pmid_dependent': res['top10_fragile1'], 'single_pmids_killing_all': res['single_pmids_killing_all_significant'],
               'fi_status_counts': dict(sorted(Counter(r['fi_status'] for r in table).items())),
               'terms_lost_without_IEP_HEP': int(sum(E['q'][i] > args.alpha for i in sig)),
               'seconds': time.perf_counter() - t0}
    (out/'summary.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    write_html(out/'report.html', summary, table)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('ledger_qc', 'inputs')}, ensure_ascii=False, indent=1))
    return 0

def write_html(path, summary, table):
    cols = ['rank', 'term', 'name', 'q', 'k', 'K', 'n_support_pmids', 'hits_without_pmid', 'single_pmid_dependent',
            'killer_pmids_direct', 'killer_pmids_via_bh', 'fragility_index', 'fi_status', 'evidence_concentration',
            'q_without_IEP_HEP', 'support_pmids']
    rows = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(r[c]))}</td>' for c in cols) + '</tr>' for r in table)
    s = {k: v for k, v in summary.items() if k not in ('ledger_qc', 'inputs')}
    page = f"""<!doctype html><html lang="en"><meta charset="utf-8"><title>pmidfrag report</title>
<style>body{{font:15px system-ui;margin:2em;max-width:1200px}}table{{border-collapse:collapse;font-size:13px}}td,th{{border:1px solid #bbb;padding:4px 6px}}th{{background:#eee;position:sticky;top:0}}</style>
<h1>pmidfrag: provenance audit of GO over-representation</h1>
<p>{html.escape(NOTE)}</p><pre>{html.escape(json.dumps(s, indent=1, ensure_ascii=False))}</pre>
<table><thead><tr>{''.join(f'<th>{c}</th>' for c in cols)}</tr></thead><tbody>{rows}</tbody></table></html>"""
    path.write_text(page, encoding='utf-8')

def main(argv=None):
    ap = argparse.ArgumentParser(prog='pmidfrag', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--version', action='version', version='pmidfrag ' + __version__)
    ap.add_argument('--gaf', required=True); ap.add_argument('--ontology', required=True)
    ap.add_argument('--background', required=True, help='one gene id per line (measured genes)')
    ap.add_argument('--query', required=True, help='one gene id per line; must be within background')
    ap.add_argument('--mapping', help='TSV accession<TAB>gene; omit if genes are the GAF object ids')
    ap.add_argument('--out', required=True)
    ap.add_argument('--db', default='UniProtKB', help='GAF column 1 database to accept (e.g. UniProtKB, SGD)')
    ap.add_argument('--taxon', default='9606'); ap.add_argument('--entity', default='protein')
    ap.add_argument('--relations', default='involved_in,acts_upstream_of,acts_upstream_of_positive_effect,acts_upstream_of_negative_effect,acts_upstream_of_or_within,acts_upstream_of_or_within_positive_effect,acts_upstream_of_or_within_negative_effect')
    ap.add_argument('--no-part-of', action='store_true', help='propagate over is_a only')
    ap.add_argument('--min-size', type=int, default=5); ap.add_argument('--max-size', type=int, default=500)
    ap.add_argument('--alpha', type=float, default=0.05)
    return run(ap.parse_args(argv))

if __name__ == '__main__':
    sys.exit(main())
