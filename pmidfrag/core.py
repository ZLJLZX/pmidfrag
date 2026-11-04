"""Task-specific provenance bounds for fixed-family overrepresentation analysis.

These are bounds over explicitly allowed annotation tables, NOT confidence
intervals, biological truth certificates, or a new FDR-control procedure.
"""
from __future__ import annotations
# Lineage: R11 robust_core (hypergeometric tail, BH, ontology projection), unchanged apart from this note.
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable
import gzip, hashlib, json, re
import numpy as np
from scipy.stats import hypergeom

UNKNOWN_REF = "urn:r11:unrecorded-reference"  # metadata flag, never a publication


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(2**20), b''): h.update(b)
    return h.hexdigest()


def bh(p: Iterable[float]) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    if p.ndim != 1 or not len(p) or not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError('BH needs a nonempty finite vector in [0,1].')
    ix = np.argsort(p, kind='stable')
    ordered = np.minimum.accumulate((p[ix] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    q = np.empty_like(p); q[ix] = np.minimum(1., ordered)
    return q


def tail(N: int, n: int, K: np.ndarray, k: np.ndarray) -> np.ndarray:
    K, k = np.asarray(K), np.asarray(k)
    if not (isinstance(N, int) and isinstance(n, int) and 0 < n <= N):
        raise ValueError('Require 0 < query size <= background size.')
    if K.shape != k.shape or K.ndim != 1 or not np.all(K == np.floor(K)) or not np.all(k == np.floor(k)):
        raise ValueError('Integer count vectors with equal shape required.')
    if np.any((K < 0) | (K > N) | (k < 0) | (k > n) | (k > K) | (K-k > N-n)):
        raise ValueError('Invalid hypergeometric counts.')
    pairs, inv = np.unique(np.column_stack([K, k]).astype(np.int64), axis=0, return_inverse=True)
    v = hypergeom.sf(pairs[:,1]-1, N, pairs[:,0], n)
    if not np.isfinite(v).all(): raise ArithmeticError('Nonfinite tail probability.')
    return v[inv]


def go_id(text: str) -> str:
    s = text.rsplit('/', 1)[-1]
    return s.replace('GO_', 'GO:')


@dataclass
class Ontology:
    terms: tuple[str, ...]
    names: dict[str, str]
    aliases: dict[str, str]
    obsolete: set[str]
    ancestors: dict[str, frozenset[str]]


def load_ontology(path: str | Path, use_part_of: bool = True) -> Ontology:
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as f: j = json.load(f)
    graph = j['graphs'][0]
    names, aliases, obsolete, bp = {}, {}, set(), set()
    alias_targets = defaultdict(set)
    for node in graph['nodes']:
        term = go_id(node['id']); meta = node.get('meta', {})
        if not re.fullmatch(r'GO:\d{7}', term): continue
        if meta.get('deprecated', False): obsolete.add(term); continue
        props = meta.get('basicPropertyValues', [])
        if not any(x['pred'].endswith('hasOBONamespace') and x['val']=='biological_process' for x in props): continue
        bp.add(term); names[term] = node.get('lbl', term)
        for x in props:
            if x['pred'].endswith('hasAlternativeId'): alias_targets[go_id(x['val'])].add(term)
    aliases = {a: next(iter(ts)) for a,ts in alias_targets.items() if len(ts)==1}
    parents = {t: set() for t in bp}
    for e in graph['edges']:
        rel = e['pred'].rsplit('/',1)[-1]
        if rel != 'is_a' and not (use_part_of and rel in ('BFO_0000050','BFO:0000050')): continue
        a, b = go_id(e['sub']), go_id(e['obj'])
        if a in bp and b in bp: parents[a].add(b)
    visiting = set()
    @lru_cache(None)
    def anc(t):
        if t in visiting: raise ValueError('Cycle in selected ontology projection.')
        visiting.add(t); out = {t}
        for p in parents[t]: out.update(anc(p))
        visiting.remove(t); return frozenset(out)
    ancestors = {t: anc(t) for t in sorted(bp)}
    root='GO:0008150'
    if root not in bp: raise ValueError('Missing biological process root.')
    terms=tuple(t for t in sorted(bp) if root in ancestors[t])
    names={t:names[t] for t in terms}; ancestors={t:ancestors[t] for t in terms}
    return Ontology(terms, names, aliases, obsolete, ancestors)


def closure(terms: Iterable[str], ancestors: dict[str,frozenset[str]]) -> frozenset[str]:
    s=set()
    for t in terms: s.update(ancestors[t])
    return frozenset(s)


def scenarios(rows, genes, ancestors, source):
    """References must be nonempty or use UNKNOWN_REF to preserve uncertainty."""
    old, exact, lower, upper = {}, {}, {}, {}
    for g in genes:
        rs=rows.get(g, [])
        if any(not r for _,r in rs): raise ValueError('Empty reference set: record missingness explicitly.')
        refs=set().union(*(r for _,r in rs)) if rs else set()
        old[g]=closure((t for t,_ in rs),ancestors)
        exact[g]=closure((t for t,r in rs if source not in r),ancestors)
        lower[g]=frozenset() if source in refs else old[g]
        upper[g]=frozenset() if refs == {source} else old[g]
        if not lower[g] <= exact[g] <= upper[g]: raise AssertionError('Annotation envelope failed.')
    return old,exact,lower,upper


def hit_counts(closed, genes, term_index):
    a=np.zeros(len(term_index),dtype=np.int64)
    for g in genes:
        ix=[term_index[t] for t in closed[g] if t in term_index]
        a[ix]+=1
    return a


def analyze(rows, genes, query, ancestors, family, source, alpha=.05, guard=1e-12):
    genes=tuple(sorted(set(genes))); query=tuple(sorted(set(query)))
    if not set(query) <= set(genes): raise ValueError('Query is not a subset of frozen background.')
    if not query: raise ValueError('Empty query: report not_tested, not a successful negative result.')
    if not family or len(family)!=len(set(family)): raise ValueError('Nonempty unique fixed family required.')
    ti={t:i for i,t in enumerate(family)}
    N,n=len(genes),len(query)
    old,exact,lower,upper=scenarios(rows,genes,ancestors,source)
    Ks=[hit_counts(s,genes,ti) for s in (old,exact,lower,upper)]
    ks=[hit_counts(s,query,ti) for s in (old,exact,lower,upper)]
    Ko,Ke,Kl,Ku=Ks; ko,ke,al,au=ks
    bl,bu=Kl-al,Ku-au
    po=tail(N,n,Ko,ko); pe=tail(N,n,Ke,ke)
    pl=tail(N,n,au+bl,au); pu=tail(N,n,al+bu,al)
    qo,qe,ql,qu=map(bh,(po,pe,pl,pu))
    qblank=bh(tail(N,n,Kl,al))
    if np.any(pl>pe+guard) or np.any(pe>pu+guard) or np.any(ql>qe+guard) or np.any(qe>qu+guard):
        raise AssertionError('Numerical envelope containment failed.')
    sig=qu < alpha-guard; no=ql > alpha+guard
    status=np.full(len(family),'unresolved',dtype='<U24')
    status[sig]='robust_significant';status[no]='robust_nonsignificant'
    eo,ee=qo<=alpha,qe<=alpha
    false=int(np.sum(sig & ~ee)+np.sum(no & ee))
    if false: raise AssertionError('A definitive call contradicts exact row-filter reanalysis.')
    interest=eo|ee
    unresolved=~(sig|no)
    guard_only=unresolved & (ql==qu)
    denom=lambda a,b: float(a/b) if b else None
    summary={
      'N':N,'n':n,'family_size':len(family),'source':source,
      'baseline_significant':int(eo.sum()),'exact_significant':int(ee.sum()),
      'unchanged_baseline_disagreements':int((eo!=ee).sum()),
      'blanket_drop_disagreements':int(((qblank<=alpha)!=ee).sum()),
      'robust_significant':int(sig.sum()),'robust_nonsignificant':int(no.sum()),
      'unresolved':int(unresolved.sum()),'unresolved_zero_width_threshold_guard':int(guard_only.sum()),
      'unresolved_nonzero_interval':int((unresolved & ~guard_only).sum()),'false_definitive':false,
      'significant_union':int(interest.sum()),'unresolved_significant_union':int((~(sig|no)&interest).sum()),
      'robust_retention_of_exact_significant':denom(sig.sum(),ee.sum()),
      'unresolved_fraction_significant_union':denom((~(sig|no)&interest).sum(),interest.sum()),
      'abstain_all_unresolved':len(family),
      'source_affected_genes':sum(any(source in r for _,r in rows.get(g,[])) for g in genes),
      'clinical_or_causal_validation':False,
    }
    arrays=dict(p_baseline=po,q_baseline=qo,p_exact=pe,q_exact=qe,p_lower=pl,p_upper=pu,
                q_lower=ql,q_upper=qu,q_blanket=qblank,K_baseline=Ko,k_baseline=ko,
                K_exact=Ke,k_exact=ke,a_lower=al,a_upper=au,b_lower=bl,b_upper=bu,interval_width_zero=(ql==qu),status=status)
    return arrays,summary
