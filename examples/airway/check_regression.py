"""Compare a pmidfrag run on the airway up-regulated list with the audited reference table.

Usage: python examples/airway/check_regression.py <pmidfrag out dir> [expected_terms_up.csv.gz] [--tamper]

Checks every term of the fixed family (K, k, q, q without IEP/HEP to 1e-12) and every significant
term (dependence flag, killer PMIDs, support count, rank), plus the internal consistency of the per-term
columns (supporting-PMID list vs count, direct + BH-only killers = all killers, FI status vs FI value).
--tamper perturbs one reference q value before comparing; the check must then FAIL (negative control).
Exit 0 = PASS, 1 = FAIL.
"""
import csv
import gzip
import json
import os
import sys

TOL = 1e-12


def main(argv):
    tamper = "--tamper" in argv
    args = [a for a in argv if a != "--tamper"]
    out = args[0]
    ref_path = args[1] if len(args) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "expected_terms_up.csv.gz")
    with gzip.open(ref_path, "rt", encoding="utf-8") as f:
        ref = {r["term"]: r for r in csv.DictReader(f)}
    if tamper:
        t = sorted(ref)[0]
        ref[t]["baseline_q"] = repr(float(ref[t]["baseline_q"]) * 0.5)
    with gzip.open(os.path.join(out, "all_terms.csv.gz"), "rt", encoding="utf-8") as f:
        got = {r["term"]: r for r in csv.DictReader(f)}
    with open(os.path.join(out, "significant_terms.csv"), encoding="utf-8") as f:
        sig = {r["term"]: r for r in csv.DictReader(f)}
    bad = []
    if set(ref) != set(got):
        bad.append("family differs: %d reference vs %d run" % (len(ref), len(got)))
    for t in sorted(set(ref) & set(got)):
        a, b = got[t], ref[t]
        if int(a["K"]) != int(b["baseline_K"]) or int(a["k"]) != int(b["baseline_k"]):
            bad.append("%s K/k" % t)
        if abs(float(a["q"]) - float(b["baseline_q"])) > TOL:
            bad.append("%s q %s vs %s" % (t, a["q"], b["baseline_q"]))
        if abs(float(a["q_without_IEP_HEP"]) - float(b["E4b1_q"])) > TOL:
            bad.append("%s q_without_IEP_HEP %s vs %s" % (t, a["q_without_IEP_HEP"], b["E4b1_q"]))
    ref_sig = {t for t, r in ref.items() if r["baseline_significant"] == "True"}
    if ref_sig != set(sig):
        bad.append("significant set differs: %d vs %d" % (len(ref_sig), len(sig)))
    for t in sorted(ref_sig & set(sig)):
        a, b = sig[t], ref[t]
        if (a["single_pmid_dependent"] == "True") != (b["fragile1"] == "True"):
            bad.append("%s dependence flag" % t)
        if set(filter(None, a["killer_pmids"].split(";"))) != set(json.loads(b["killer_pmids"])):
            bad.append("%s killer PMIDs" % t)
        if int(a["n_support_pmids"]) != int(b["support_pmids"]):
            bad.append("%s support" % t)
        if len(set(filter(None, a["support_pmids"].split(";")))) != int(a["n_support_pmids"]):
            bad.append("%s support list vs count" % t)
        split = [set(filter(None, a[c].split(";"))) for c in ("killer_pmids", "killer_pmids_direct", "killer_pmids_via_bh")]
        if split[0] != split[1] | split[2] or split[1] & split[2]:
            bad.append("%s killer split" % t)
        fi, st = a["fragility_index"], a["fi_status"]
        if not ((fi in ("1", "2") and st == "exact") or (fi in ("3", "4", "5") and st == "upper_bound")
                or (fi == "" and st.startswith("not_abolished"))):
            bad.append("%s FI %r with status %r" % (t, fi, st))
        if int(a["rank"]) != int(b["significant_rank"]):
            bad.append("%s rank" % t)
    print("[check] family terms compared: %d; significant terms compared: %d" % (len(set(ref) & set(got)), len(ref_sig & set(sig))))
    for x in bad[:20]:
        print("[check] MISMATCH", x)
    print("[check] %s (%d mismatches)%s" % ("PASS" if not bad else "FAIL", len(bad), " [tampered reference]" if tamper else ""))
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
