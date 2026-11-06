"""GO biological-process ontology loading, propagation, hypergeometric tail and BH adjustment."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable
import gzip, hashlib, json, re
import numpy as np
from scipy.stats import hypergeom


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


def hit_counts(closed, genes, term_index):
    a=np.zeros(len(term_index),dtype=np.int64)
    for g in genes:
        ix=[term_index[t] for t in closed[g] if t in term_index]
        a[ix]+=1
    return a
