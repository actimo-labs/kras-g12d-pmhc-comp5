#!/usr/bin/env python
"""Pick 12 pilot designs spanning the iptm_margin_wt range, for the sample test.

Spanning the range matters: if the test only sampled designs clustered near
zero, a rank correlation computed on them would be dominated by noise regardless
of whether averaging helped. Spread the sample across the observable range so
the correlation has something to resolve.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
rows = list(csv.DictReader(open(ROOT / "pilot/boltz/main_three_state.csv")))


def g(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return 0.0


rows.sort(key=lambda r: g(r, "iptm_margin_wt_peptide"))
n = len(rows)
idx = sorted({0, 1, 2, n // 4, n // 4 + 1, n // 2, n // 2 + 1,
              3 * n // 4, 3 * n // 4 + 1, n - 3, n - 2, n - 1})
print(",".join(rows[i]["candidate_id"] for i in idx if 0 <= i < n))
