import yaml, sys

CFG = "/mnt/user-data/uploads/configs/challenge.yaml"
cfg = yaml.safe_load(open(CFG))

VALID = set("ACDEFGHIKLMNPQRSTVWY")
ok = True

def check(name, seq, claimed=None):
    global ok
    bad = sorted(set(seq) - VALID)
    lenok = (claimed is None) or (len(seq) == claimed)
    if bad or not lenok:
        ok = False
    print(f"  {name:16s} len={len(seq):4d}"
          + (f"  claimed={claimed:4d}  match={lenok}" if claimed else "")
          + f"  invalid={bad or 'none'}")

t = cfg["target"]
print("TARGET CHAINS (from challenge.yaml)")
check("heavy_chain", t["heavy_chain"]["seq"], t["heavy_chain"]["length"])
check("b2m",         t["b2m"]["seq"],         t["b2m"]["length"])

print("\nPEPTIDES")
g12d = cfg["states"]["g12d"]["peptide"]
wt   = cfg["states"]["wt"]["peptide"]
check("g12d", g12d, 10)
check("wt",   wt,   10)

# UniProt range consistency
hr, br = t["heavy_chain"]["range"], t["b2m"]["range"]
pr = cfg["states"]["g12d"]["peptide_range"]
for nm, rng, seq in [("heavy_chain", hr, t["heavy_chain"]["seq"]),
                     ("b2m", br, t["b2m"]["seq"]),
                     ("peptide", pr, g12d)]:
    span = rng[1] - rng[0] + 1
    flag = "OK" if span == len(seq) else "MISMATCH"
    print(f"  range {nm:12s} {rng[0]}-{rng[1]} = {span} residues vs seq {len(seq)}  [{flag}]")
    if span != len(seq):
        ok = False

print("\nMUTATION CHECK")
diffs = [(i+1, a, b) for i, (a, b) in enumerate(zip(g12d, wt)) if a != b]
print(f"  differences (1-indexed): {diffs}")
exp = [(6, 'D', 'G')]
print(f"  expected: {exp}  -> {'OK' if diffs == exp else 'MISMATCH'}")
if diffs != exp:
    ok = False

print("\nHOTSPOT SANITY (peptide chain C)")
hs = cfg["hotspots"]
for p in hs["primary"]:
    print(f"  p{p:<2d} = {g12d[p-1]}  (primary hotspot)")
for p in hs["tracked_not_targeted"]:
    print(f"  p{p:<2d} = {g12d[p-1]}  (tracked, buried in bulge_down)")
for p in hs["anchors_excluded"]:
    print(f"  p{p:<2d} = {g12d[p-1]}  (anchor, excluded)")
covered = sorted(hs["primary"] + hs["tracked_not_targeted"] + hs["anchors_excluded"])
missing = [i for i in range(1, len(g12d)+1) if i not in covered]
print(f"  positions not classified: {missing}")

print("\nCYSTEINES (disulfides to preserve in prep)")
print(f"  heavy_chain Cys: {[i+1 for i,c in enumerate(t['heavy_chain']['seq']) if c=='C']}")
print(f"  b2m Cys:         {[i+1 for i,c in enumerate(t['b2m']['seq']) if c=='C']}")
print(f"  peptide Cys:     {[i+1 for i,c in enumerate(g12d) if c=='C'] or 'none'}")

print("\nSUBMISSION GUARDS")
ch = cfg["challenge"]
print(f"  track={ch['track']}  max_designs={ch['max_designs']}  len=[{ch['length_min']},{ch['length_max']}]")
print(f"  submission.columns_verified = {cfg['submission']['columns_verified']}")
print(f"  auto_submit = {cfg['submission']['auto_submit']}")
prov = [k for k, v in cfg["thresholds"].items() if isinstance(v, dict) and v.get("provisional")]
print(f"  provisional thresholds still to calibrate: {prov}")

print("\n" + ("ALL SEQUENCE CHECKS PASSED" if ok else "*** CHECKS FAILED — DO NOT PROCEED ***"))

# emit FASTA for the three assay states
out = "/tmp/claude-0/-home-claude/a45be8bd-58e1-50e8-9eb3-6cc2a283c48a/scratchpad/target_states.fasta"
hc, b2m = t["heavy_chain"]["seq"], t["b2m"]["seq"]
with open(out, "w") as f:
    f.write(f">g12d_A_heavy|HLA-A*11:01|Q5S3G3_25-305\n{hc}\n")
    f.write(f">g12d_B_b2m|P61769_21-119\n{b2m}\n")
    f.write(f">g12d_C_peptide|KRAS_G12D|P01116_7-16\n{g12d}\n")
    f.write(f">wt_A_heavy|HLA-A*11:01|Q5S3G3_25-305\n{hc}\n")
    f.write(f">wt_B_b2m|P61769_21-119\n{b2m}\n")
    f.write(f">wt_C_peptide|KRAS_WT|P01116_7-16\n{wt}\n")
    f.write(f">empty_A_heavy|HLA-A*11:01|Q5S3G3_25-305\n{hc}\n")
    f.write(f">empty_B_b2m|P61769_21-119\n{b2m}\n")
print(f"\nwrote {out}")
sys.exit(0 if ok else 1)