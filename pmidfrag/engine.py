"""E4 PLAN + amendments 1/2, implemented without species-specific data selection.

Accepted rows, background, query, ancestor table and fixed family are explicit
inputs. A PMID deletion removes rows from the ENTIRE background. Non-query-only
sources cannot increase any overrepresentation p value, or any fixed-family BH
q value, and thus cannot kill a baseline-significant term.
"""
from __future__ import annotations
from collections import Counter, defaultdict
import re
from typing import Callable
import numpy as np
from scipy.stats import hypergeom
from .core import bh, tail
from .ledger import EXPRESSION_CODES

PMID = re.compile(r'PMID:[0-9]+', re.ASCII)
PRIMARY_SEED = 20261102
PRIMARY_REPLICATES = 200
HEP_SOURCE = 'PMID:22926577'
HEP_TERM = 'GO:0021762'


def pmids(references):
    return frozenset(x for x in references if isinstance(x, str) and PMID.fullmatch(x))


def primary_control_plan(records, background, repetitions=PRIMARY_REPLICATES,
                         seed=PRIMARY_SEED):
    """A3: one RNG across replicate-major, lexicographic-gene traversal.

    Within a gene, rows retain accepted-ledger order. A shortage removes ALL
    eligible nonexpression rows, records the deficit, and does not consume RNG
    draws. At equality, choice(..., replace=False) is called as for sufficient
    pools. IDs are sorted only for serialization, AFTER sampling.
    """
    if not isinstance(repetitions, int) or repetitions <= 0:
        raise ValueError('positive repetitions required')
    bg = set(background)
    groups = {g: [] for g in sorted(bg)}
    ids = [r['row_id'] for r in records]
    if len(ids) != len(set(ids)):
        raise ValueError('Accepted row IDs must be unique, including duplicate records')
    for r in records:
        if r['gene'] not in bg:
            raise ValueError('Accepted record lies outside frozen background')
        groups[r['gene']].append(r)
    specs, pools, target_ids = [], {}, []
    for gene in sorted(bg):
        target = [r['row_id'] for r in groups[gene] if r['evidence'] in EXPRESSION_CODES]
        pool = [r['row_id'] for r in groups[gene] if r['evidence'] not in EXPRESSION_CODES]
        k = len(target)
        if k:
            pools[gene] = (pool, k)
            target_ids.extend(target)
        specs.append({'gene': gene, 'target_rows': k, 'control_pool': len(pool),
                      'removed_per_replicate': min(k, len(pool)),
                      'shortfall_rows': max(0, k-len(pool)),
                      'status': ('no_target' if not k else
                                 ('shortfall_delete_all' if len(pool)<k else 'matched'))})
    rng = np.random.default_rng(seed)
    draws = []
    for rep in range(repetitions):
        removed = []
        for gene in sorted(pools):
            pool, k = pools[gene]
            if len(pool) < k:
                removed.extend(pool)
            else:
                selected = rng.choice(len(pool), size=k, replace=False)
                removed.extend(pool[int(i)] for i in selected)
        if len(removed) != len(set(removed)):
            raise AssertionError('Repeated row within one deletion draw')
        draws.append(sorted(removed))
    total = len(target_ids)
    actual = sum(s['removed_per_replicate'] for s in specs)
    deficit = sum(s['shortfall_rows'] for s in specs)
    assert actual + deficit == total
    return {'kind': 'E4_A3_PRIMARY', 'replicates': repetitions, 'seed': seed,
            'rng': 'numpy.random.default_rng', 'bit_generator': type(rng.bit_generator).__name__,
            'numpy_version': np.__version__,
            'traversal': 'replicate outer; gene_id lexicographic inner; ledger-order pool',
            'shortfall_policy': 'delete_all_available_and_record_deficit; no_rng_draw',
            'equal_pool_policy': 'choice_without_replacement',
            'per_gene': specs, 'expression_row_ids': sorted(target_ids),
            'control_row_ids': draws, 'target_total': total,
            'removed_rows_per_replicate': actual, 'shortfall_rows_per_replicate': deficit,
            'genes_with_shortfall': sum(s['shortfall_rows'] > 0 for s in specs),
            'target_row_coverage': actual/total if total else None,
            'status': 'ready' if total else 'no_target_expression_rows',
            'independent_biological_replicates': 0}


def validate_primary_plan(plan, records, background):
    if plan.get('kind') != 'E4_A3_PRIMARY' or plan.get('seed') != PRIMARY_SEED or plan.get('replicates') != PRIMARY_REPLICATES:
        raise ValueError('A3 specification drift')
    expected = primary_control_plan(records, background)
    if plan != expected:
        raise ValueError('A3 row IDs, RNG environment, counts or shortfalls differ from lock')
    return True


class FragilityEngine:
    """Exact sparse update of known row deletion; full fixed-family BH each time.

    No source is selected on enrichment outcomes. Candidates include ALL valid
    query-row PMIDs even when their deletion is certified to change no closure.
    """
    def __init__(self, records, background, ancestors, family):
        self.bg = tuple(sorted(set(background)))
        self.bgset = set(self.bg)
        self.family = tuple(family)
        if not self.bg or not self.family or len(self.family) != len(set(self.family)):
            raise ValueError('Nonempty background and unique fixed family required')
        self.ti = {t:i for i,t in enumerate(self.family)}
        self.anc = ancestors
        self.rows = {g: [] for g in self.bg}
        self.sources = defaultdict(set)
        self.records = list(records)
        for r in self.records:
            if r['gene'] not in self.bgset or r['term'] not in ancestors:
                raise ValueError('Record outside background or missing ancestor definition')
            if not r['references']:
                raise ValueError('Missing references need an explicit metadata flag')
            self.rows[r['gene']].append(r)
            for s in pmids(r['references']):
                self.sources[s].add(r['gene'])
        self.term_ix = {
            t: frozenset(self.ti[a] for a in ancestors[t] if a in self.ti)
            for t in {r['term'] for r in self.records}}
        self.closed = {g: self._close(self.rows[g]) for g in self.bg}
        self.K = self._count(self.bg)
        self.delta_cache = {}
        self.tail_cache = {}

    def _close(self, rows):
        out = set()
        for r in rows:
            out.update(self.term_ix[r['term']])
        return frozenset(out)

    def _count(self, genes):
        v = np.zeros(len(self.family), dtype=np.int64)
        for g in genes:
            if self.closed[g]:
                v[list(self.closed[g])] += 1
        return v

    def query_set(self, query):
        q = set(query)
        if not q <= self.bgset:
            raise ValueError('Query outside frozen background')
        return tuple(sorted(q))

    def _tail(self, n, K, k):
        """R10/E2-style exact memoization; no statistical approximation."""
        N=len(self.bg);K=np.asarray(K);k=np.asarray(k)
        if np.any((K<0)|(K>N)|(k<0)|(k>n)|(k>K)|(K-k>N-n)):
            raise ValueError('Invalid hypergeometric counts')
        result=np.ones(len(K),dtype=float)
        ids=np.flatnonzero(k>np.maximum(0,n+K-N))
        if not len(ids):return result
        pairs,inv=np.unique(np.column_stack([K[ids],k[ids]]),axis=0,return_inverse=True)
        missing=[(n,int(a),int(b)) for a,b in pairs if (n,int(a),int(b)) not in self.tail_cache]
        if missing:
            ar=np.asarray(missing)
            values=hypergeom.sf(ar[:,2]-1,N,ar[:,1],ar[:,0])
            if not np.isfinite(values).all():raise ArithmeticError('Nonfinite tail')
            self.tail_cache.update(zip(missing,map(float,values)))
        result[ids]=np.asarray([self.tail_cache[(n,int(a),int(b))] for a,b in pairs])[inv]
        return result

    def baseline(self, query):
        q = self.query_set(query)
        if not q:
            raise ValueError('Empty query is not tested')
        k = self._count(q)
        p = self._tail(len(q), self.K, k)
        return {'K': self.K.copy(), 'k': k, 'p': p, 'q': bh(p)}

    def candidates(self, query):
        return sorted({s for g in self.query_set(query) for r in self.rows[g]
                       for s in pmids(r['references'])})

    def source_delta(self, source):
        if source in self.delta_cache:
            return self.delta_cache[source]
        lost = {}
        for g in sorted(self.sources.get(source, ())):
            kept = self._close(r for r in self.rows[g] if source not in r['references'])
            gone = self.closed[g] - kept
            if gone:
                lost[g] = np.asarray(sorted(gone), dtype=np.int64)
        if lost:
            ix, counts = np.unique(np.concatenate(list(lost.values())), return_counts=True)
        else:
            ix, counts = np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
        answer = (lost, ix, counts)
        self.delta_cache[source] = answer
        return answer

    def remove_source(self, query, source, baseline=None):
        q = self.query_set(query)
        if not q:
            raise ValueError('Empty query is not tested')
        base = self.baseline(q) if baseline is None else baseline
        lost, ix, diff = self.source_delta(source)
        K, k, p = base['K'].copy(), base['k'].copy(), base['p'].copy()
        if len(ix):
            K[ix] -= diff
            for g in q:
                if g in lost:
                    k[lost[g]] -= 1
            # All query changes are among indices where the background changed.
            p[ix] = self._tail(len(q), K[ix], k[ix])
        return {'K': K, 'k': k, 'p': p, 'q': bh(p)}

    def support(self, query):
        m = len(self.family)
        supporting_pmids = [set() for _ in range(m)]
        no_pmid = np.zeros(m, dtype=np.int64)
        expression_only = np.zeros(m, dtype=np.int64)
        specific_hep_only = np.zeros(m, dtype=np.int64)
        for g in self.query_set(query):
            any_pmid, non_expression, other_than_specific_hep = set(), set(), set()
            for r in self.rows[g]:
                ix = self.term_ix[r['term']]
                ids = pmids(r['references'])
                if ids:
                    any_pmid.update(ix)
                    for j in ix:
                        supporting_pmids[j].update(ids)
                if r['evidence'] not in EXPRESSION_CODES:
                    non_expression.update(ix)
                if not (r['evidence'] == 'HEP' and HEP_SOURCE in r['references']):
                    other_than_specific_hep.update(ix)
            if self.closed[g]:
                no_pmid[list(self.closed[g] - any_pmid)] += 1
                expression_only[list(self.closed[g] - non_expression)] += 1
                specific_hep_only[list(self.closed[g] - other_than_specific_hep)] += 1
        return {'support_pmid_ids': [sorted(s) for s in supporting_pmids],
                'support_pmids': np.asarray([len(s) for s in supporting_pmids], dtype=np.int64),
                'only_no_pmid_hit_genes': no_pmid,
                'only_IEP_HEP_hit_genes': expression_only,
                'only_target_HEP_hit_genes': specific_hep_only}

    def evaluate(self, query, alpha=.05, source_sink: Callable | None = None):
        q = self.query_set(query)
        if not q:
            return {'status': 'not_tested_empty_query', 'N': len(self.bg), 'n': 0,
                    'candidate_pmids': [], 'baseline': None}
        base = self.baseline(q)
        sig = np.flatnonzero(base['q'] <= alpha)
        candidates = self.candidates(q)
        killers = [[] for _ in self.family]
        source_counts, identity_count = {}, 0
        # Record every candidate, including closures with no change.
        qset=set(q)
        for source in candidates:
            delta=self.source_delta(source)
            # With no requested score stream, use the same safe screening as E2:
            # no query family hit can be lost => p/q can only stay/decrease.
            # Still record the candidate and its exact zero-killer outcome.
            if source_sink is None and (not len(sig) or not (qset & delta[0].keys())):
                source_counts[source]=0
                identity_count += int(not len(delta[1]))
                continue
            exact = self.remove_source(q, source, base)
            lost_terms = sig[exact['q'][sig] > alpha]
            for j in lost_terms:
                killers[int(j)].append(source)
            source_counts[source] = len(lost_terms)
            identity_count += int(not len(self.source_delta(source)[1]))
            if source_sink is not None:
                source_sink(source, exact)
        order = sorted(map(int,sig), key=lambda j: (float(base['q'][j]), -int(base['k'][j]), self.family[j]))
        top = order[:10]
        fragile = np.asarray([bool(s) for s in killers])
        all_killers = sorted(s for s,c in source_counts.items() if len(sig) and c == len(sig))
        support = self.support(q)
        no_pmid_rows = sum(not pmids(r['references']) for g in q for r in self.rows[g])
        return {'status': 'computed' if len(sig) else 'not_tested_no_baseline_significant',
                'N': len(self.bg), 'n': len(q), 'family_size': len(self.family),
                'baseline': base, 'baseline_significant_indices': order, 'top_indices': top,
                'candidate_pmids': candidates, 'candidate_count': len(candidates),
                'candidate_identity_count': identity_count,
                'query_accepted_rows_without_pmid': no_pmid_rows,
                'killers': killers, 'fragile1': fragile, 'support': support,
                'baseline_significant': len(sig), 'fragile1_count': int(fragile.sum()),
                'fragile1_fraction': float(fragile.sum()/len(sig)) if len(sig) else None,
                'top_term': self.family[order[0]] if order else None,
                'top_term_fragile1': bool(fragile[order[0]]) if order else None,
                'top10_count': len(top), 'top10_fragile1': int(fragile[top].sum()),
                'single_pmids_killing_all_significant': all_killers,
                'source_killed_counts': source_counts}


def iep_composition(records, names):
    rows = [r for r in records if r['evidence']=='IEP']
    counts = Counter(r['term'] for r in rows)
    missing = sorted(t for t in counts if t not in names)
    if missing:
        raise ValueError('Direct term names missing for IEP composition: '+repr(missing[:5]))
    n_response = sum(names[r['term']].startswith(('response to','cellular response to')) for r in rows)
    return {'scope':'accepted IEP rows mapped to frozen measured background; direct terms, not ancestors',
            'IEP_rows':len(rows), 'response_prefix_rows':n_response,
            'response_prefix_fraction':n_response/len(rows) if rows else None,
            'top20_direct_terms':[{'term':t,'term_name':names[t],'rows':c}
                                  for t,c in sorted(counts.items(),key=lambda kv:(-kv[1],kv[0]))[:20]],
            'descriptive_only':True}


def delta_for_row_ids(engine, remove_ids):
    """Sparse deletion for a locked row-ID set; all background rows are eligible."""
    ids = set(remove_ids)
    all_ids = {r['row_id'] for r in engine.records}
    if not ids <= all_ids:
        raise ValueError('Deletion refers to unknown accepted row IDs')
    affected = {r['gene'] for r in engine.records if r['row_id'] in ids}
    lost = {}
    for g in sorted(affected):
        kept = engine._close(r for r in engine.rows[g] if r['row_id'] not in ids)
        gone = engine.closed[g] - kept
        if gone:
            lost[g] = np.asarray(sorted(gone), dtype=np.int64)
    if lost:
        ix, count = np.unique(np.concatenate(list(lost.values())), return_counts=True)
    else:
        ix, count = np.empty(0,dtype=np.int64), np.empty(0,dtype=np.int64)
    return lost, ix, count


def score_delta(engine, query, baseline, delta):
    qgenes = engine.query_set(query)
    lost, ix, count = delta
    K,k,p = baseline['K'].copy(), baseline['k'].copy(), baseline['p'].copy()
    if len(ix):
        K[ix] -= count
        for g in qgenes:
            if g in lost:
                k[lost[g]] -= 1
        p[ix] = engine._tail(len(qgenes),K[ix],k[ix])
    return {'K':K,'k':k,'p':p,'q':bh(p)}
