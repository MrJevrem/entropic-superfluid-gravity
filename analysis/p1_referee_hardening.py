#!/usr/bin/env python
"""
Paper I referee-hardening:
 (1) numeric checks specific to Paper I (acoustic-Planck scales; M_Pl,ac notation trap)
 (2) PORTFOLIO LINT: scan all five drafts for claims deprecated by the v0.2 hardening
     passes (stale comparators, the Omega_L=0.70 artifacts, bar overclaim, protocol
     overclaim, stale counts). Re-runnable before submission.
"""
import os, sys, re, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2_guard import DERIVED, ROOT

print("=== (1) acoustic-Planck scales (m = 8.5 eV, galactic c_s) ===")
EV = 1.602176634e-19; HBAR = 1.054571817e-34; C = 2.99792458e8
m_kg = 8.5 * EV / C**2
for vf in (150e3, 220e3):
    cs = vf / np.sqrt(2)
    xi = HBAR / (m_kg * cs)
    print(f"v_f={vf/1e3:.0f} km/s: xi = {xi*1e3:.3f} mm, l_Pl,ac = sqrt(12pi)*xi = {np.sqrt(12*np.pi)*xi*1e3:.2f} mm")
print(f"M_Pl,ac = m/sqrt(12pi) = 8.5/{np.sqrt(12*np.pi):.3f} = {8.5/np.sqrt(12*np.pi):.2f} eV"
      f"   <- draft writes '0.163 m' (reads as METERS; must be '~1.4 eV' or '0.163m')")

print("\n=== (2) portfolio lint: claims deprecated by v0.2 hardening ===")
DOCS = ["docs/papers/PAPER_I_v2.md", "docs/papers/PAPER_II_v2.md", "docs/papers/PAPER_III_v2.md",
        "docs/papers/PAPER_IV_v2.md", "docs/papers/PAPER_V_v2.md", "docs/papers/LETTER_v2.md",
        "docs/outreach/PITCH.md", "docs/outreach/FAQ.md", "docs/outreach/TWEETS.md"]
# pattern -> why it is deprecated (changelog mentions are fine; flag for human judgment)
RULES = [
    (r"69\.6",                 "H0=69.6 was the Omega_L=0.70 artifact (Letter v0.2: 70.3-74.8 degeneracy)"),
    (r"0\.13\s*%",             "1/6-proximity 0.13% was at Omega_L=0.70; Planck value is ~1%"),
    (r"\+4\.4\s*%(?![^.]*11)", "sole-comparator +4.4%; v0.2 standard is '4-11% below' both comparators"),
    (r"bars (stay|are) fast|fast bars|bars fast", "bar is Mach 1.4 supersonic: reduced-not-zero drag (P3 v0.2)"),
    (r"cryptograph",           "protocol overclaim; v0.2 standard: commit-stamped + public deposit"),
    (r"git-verified",          "same protocol overclaim"),
    (r"occupancy [^L]*46[06]", "occupancy now the Liouville pair ~500 (503 cosmic / 467 halo)"),
    (r"three methods and two arrays", "P2 v0.2 demoted the three-way invariance (budget-convention coincidence)"),
    (r"52 commits",            "stale commit count; use generic wording"),
    (r"1\.128[^ ]* ?m s|value 1\.128", "1.128 as THE comparator; canonical 1.20 must be primary (Letter v0.2)"),
    # --- appended by the Paper IV pass ---
    (r"new exact (NLS )?solution", "novelty overclaim: 1D TG antecedent (Kolomeisky+ 2000); claim 3D realization/stability instead (P4 v0.2)"),
    (r"attractive (light )?sheet", "pinning sign error: dark solitons pin on BLUE-detuned/repulsive barriers (P4 v0.2)"),
    # --- appended by the figure pass ---
    (r"slope √\(2/3\) vs 1/√2|steeper.*√\(2/3\)", "convention mix: at matched xi=hbar/mc_s the cubic slope is 1; quintic is SHALLOWER/broader (figure pass)"),
    # --- appended by the authorship/formatting pass ---
    (r"arXiv (posting|deposit)", "self-referential arXiv mention; use 'public deposit/posting' (author pass)"),
    (r"[Nn]ot peer-reviewed", "working-draft disclaimer removed for public posting (author pass)"),
    (r"Bullet Cluster to Gaia DR4", "Paper III retitled — no future-data claim in a title (author pass)"),
    (r"References \(draft\)", "section renamed 'References' for public posting (author pass)"),
    # --- appended by the T3/Paper V sync pass ---
    (r"500–800 km|500-800 km|σ_crit ≈ 500", "narrative σ_crit band superseded by T3/Paper V: kernel window 250–1000 (ρ-spread 476–756 at m_fid); break measured, location systematics-limited"),
    (r"groups straddle", "superseded by T3: X-ray groups uniformly elevated; transition sits between field-ETG and group scales"),
    # --- appended by the authorship-decision pass (2026-07-25) ---
    (r"and Claude Fable 5", "AI co-author line removed per arXiv/APS policy; credit = byline subline + AI contribution statement (user decision 2026-07-25)"),
]
hits = []
for doc in DOCS:
    path = os.path.join(ROOT, doc)
    for ln, line in enumerate(open(path, encoding="utf-8"), 1):
        for pat, why in RULES:
            if re.search(pat, line, re.I):
                ctx = line.strip()[:90]
                hits.append((doc, ln, pat, why, ctx))
if hits:
    cur = None
    for doc, ln, pat, why, ctx in hits:
        if doc != cur:
            print(f"\n[{doc}]"); cur = doc
        print(f"  L{ln:4d}  {why}\n        > {ctx}")
else:
    print("clean")
print(f"\n{len(hits)} flagged lines (changelog/self-referential mentions are acceptable; judge each)")
json.dump([dict(doc=d, line=l, why=w) for d, l, _, w, _ in hits],
          open(os.path.join(DERIVED, "p1_portfolio_lint.json"), "w"), indent=1)
print("saved derived/p1_portfolio_lint.json")