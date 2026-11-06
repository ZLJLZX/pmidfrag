"""Lossless GAF row ledger: every row is kept with its inclusion or exclusion reason.

Every nonblank, nonheader GAF row has a unique physical-line identity and a
terminal inclusion/exclusion reason. Duplicate records are NOT deduplicated.
Column 6 identifiers are retained as aliases on ONE reference-bearing row;
there is no claim that each pipe-separated identifier is an independent study.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, re
from collections import Counter
from pathlib import Path
from typing import Iterable

UNKNOWN_REF = 'urn:pmidfrag:unrecorded-reference'  # metadata flag, never a publication
EXPRESSION_CODES = frozenset({'IEP', 'HEP'})
KNOWN_EVIDENCE_CODES = frozenset(('EXP IDA IPI IMP IGI IEP HTP HDA HMP HGI HEP IBA IBD IKR IRD ISS ISO ISA ISM IGC RCA TAS NAS IC ND IEA').split())


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def classify(z: list[str], p2g: dict[str, str], ontology, relations: set[str], db: str = 'UniProtKB', taxon: str = '9606', entity: str = 'protein'):
    """Return (reason, normalized row); exactly one terminal reason per GAF row."""
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


def filter_records(records, source=None, exclude_codes=(), remove_ids=()):
    codes, removed = set(exclude_codes), set(remove_ids)
    return [r for r in records if r['evidence'] not in codes
            and (source is None or source not in r['references'])
            and r['row_id'] not in removed]
