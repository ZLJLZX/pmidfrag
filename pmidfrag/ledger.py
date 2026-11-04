"""Add a lossless row ledger without changing the R11 inclusion or scoring policy.

Every nonblank, nonheader GAF row has a unique physical-line identity and a
terminal inclusion/exclusion reason. Duplicate records are NOT deduplicated.
Column 6 identifiers are retained as aliases on ONE reference-bearing row;
there is no claim that each pipe-separated identifier is an independent study.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

UNKNOWN_REF = 'urn:r11:unrecorded-reference'
EXPRESSION_CODES = frozenset({'IEP', 'HEP'})
KNOWN_EVIDENCE_CODES = frozenset(('EXP IDA IPI IMP IGI IEP HTP HDA HMP HGI HEP IBA IBD IKR IRD ISS ISO ISA ISM IGC RCA TAS NAS IC ND IEA').split())


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def classify(z: list[str], p2g: dict[str, str], ontology, relations: set[str], db: str = 'UniProtKB', taxon: str = '9606', entity: str = 'protein'):
    """Mirror frozen R11 parse_gaf; do not silently introduce new exclusions."""
    if len(z) != 17: return 'malformed_column_count', None
    if z[0] != db: return 'other_db', None
    if z[8] != 'P': return 'not_BP', None
    if z[12].split('|')[0].split(':')[-1] != taxon: return 'other_taxon', None
    if entity and z[11] != entity: return 'other_entity', None
    q = set(z[3].split('|'))
    if 'NOT' in q: return 'NOT', None
    if len(q) != 1 or not q <= relations: return 'other_relation', None
    if z[6] == 'ND': return 'ND', None
    if z[15]: return 'annotation_extension', None
    protein = re.sub(r'-\d+$', '', z[1].removeprefix(db + ':'))
    if protein not in p2g: return 'unmapped_or_outside_background', None
    term = z[4]
    if term in ontology.obsolete: return 'obsolete_term', None
    if term not in ontology.ancestors: term = ontology.aliases.get(term, term)
    if term not in ontology.ancestors: return 'unresolved_or_nonBP_term', None
    refs = sorted({r.strip() for r in z[5].split('|') if r.strip()})
    return 'accepted', {
        'gene': p2g[protein], 'protein': protein, 'term': term,
        'references': refs or [UNKNOWN_REF], 'evidence': z[6],
        'reference_missing': not bool(refs),
    }


def build_ledger(gaf: str | Path, p2g: dict[str, str], ontology,
                 relations: Iterable[str], out: str | Path, db: str = 'UniProtKB', taxon: str = '9606', entity: str = 'protein'):
    out = Path(out)
    if out.exists(): raise FileExistsError(f'Refuse to overwrite ledger: {out}')
    out.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(gaf).endswith('.gz') else open
    accepted = []; reasons = Counter(); evidence = Counter(); header = []
    nonblank = 0
    with opener(gaf, 'rb') as f, gzip.open(out, 'wt', encoding='utf-8', newline='\n') as w:
        for line_number, raw in enumerate(f, 1):
            text = raw.decode('utf-8').rstrip('\r\n')
            if text.startswith('!'):
                header.append({'line': line_number, 'text': text}); continue
            if not text.strip(): continue
            nonblank += 1
            fields = text.split('\t')
            reason, norm = classify(fields, p2g, ontology, set(relations), db, taxon, entity)
            digest = hashlib.sha256(raw).hexdigest()
            row_id = f'{line_number}:{digest}'
            row = {'row_id': row_id, 'line_number': line_number,
                   'line_sha256_with_terminator': digest,
                   'fields': fields, 'status': reason, 'normalized': norm}
            w.write(json.dumps(row, ensure_ascii=False, separators=(',', ':'))+'\n')
            reasons[reason] += 1
            if norm is not None:
                item = {'row_id': row_id, 'line_number': line_number, **norm}
                accepted.append(item); evidence[norm['evidence']] += 1
    if sum(reasons.values()) != nonblank: raise AssertionError('Row ledger does not balance')
    summary = {'input_sha256': sha256(gaf), 'input_rows': nonblank,
               'reason_counts': dict(sorted(reasons.items())),
               'accepted_evidence_counts': dict(sorted(evidence.items())),
               'accepted_rows': len(accepted), 'headers': header,
               'ledger_sha256': sha256(out)}
    return accepted, summary


def project(records):
    out = defaultdict(list)
    for r in records:
        out[r['gene']].append((r['term'], frozenset(r['references'])))
    return dict(out)


def assert_projection(records, expected):
    got = project(records)
    want = {g: [(t, frozenset(refs)) for t, refs in rs] for g, rs in expected.items() if rs}
    if got != want:
        bad = [g for g in sorted(set(got) | set(want)) if got.get(g) != want.get(g)]
        raise AssertionError(f'Projection differs from frozen R11 for {len(bad)} genes; first={bad[:5]}')
    return {'genes_with_rows': len(got), 'rows': sum(map(len,got.values())),
            'ordered_rows_and_multiplicity_identical': True}


def filter_records(records, source=None, exclude_codes=(), remove_ids=()):
    codes, removed = set(exclude_codes), set(remove_ids)
    return [r for r in records if r['evidence'] not in codes
            and (source is None or source not in r['references'])
            and r['row_id'] not in removed]


def matched_control_plan(records, repetitions=20, salt='R13-nonexpression-control-v1|'):
    """Same-gene, same-row-count deletion. Unmatchable targets are explicit.

    Hash-ranked rows give deterministic pseudo-random controls, not biological
    replicates. All deletions must be locked before differential expression.
    """
    if not isinstance(repetitions, int) or repetitions <= 0: raise ValueError('positive repetitions required')
    groups = defaultdict(list)
    for r in records: groups[r['gene']].append(r)
    matched, unmatched, manifest, control_sets = [], [], [], [[] for _ in range(repetitions)]
    for gene in sorted(groups):
        target = [r for r in groups[gene] if r['evidence'] in EXPRESSION_CODES]
        if not target: continue
        pool = [r for r in groups[gene] if r['evidence'] not in EXPRESSION_CODES]
        n = len(target)
        if len(pool) < n:
            unmatched.extend(r['row_id'] for r in target)
            manifest.append({'gene': gene, 'target_rows': n, 'control_pool': len(pool), 'status': 'unmatchable'})
            continue
        matched.extend(r['row_id'] for r in target)
        manifest.append({'gene': gene, 'target_rows': n, 'control_pool': len(pool), 'status': 'matched'})
        for rep in range(repetitions):
            ordered = sorted(pool, key=lambda r: hashlib.sha256(f'{salt}{rep}|{gene}|{r["row_id"]}'.encode()).hexdigest())
            control_sets[rep].extend(r['row_id'] for r in ordered[:n])
    for ids in control_sets:
        if len(ids) != len(matched) or len(ids) != len(set(ids)):
            raise AssertionError('Matched-control identity/count violation')
    target_total = len(matched) + len(unmatched)
    return {'replicates': repetitions, 'salt': salt, 'per_gene': manifest,
            'matched_expression_row_ids': sorted(matched),
            'unmatched_expression_row_ids': sorted(unmatched),
            'control_row_ids': [sorted(ids) for ids in control_sets],
            'target_row_coverage': len(matched)/target_total if target_total else None,
            'target_total': target_total, 'matched_target_count': len(matched),
            'status': 'no_target_expression_rows' if not target_total else ('no_matched_genes' if not matched else 'ready'),
            'independent_biological_replicates': 0}
