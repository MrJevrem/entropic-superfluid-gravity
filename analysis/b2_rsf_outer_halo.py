#!/usr/bin/env python
"""
D24 — the R_SF outer-halo test (BK-comparison open item 3).
Question: our MW superfluid extends to R_SF ~ 185-310 kpc (RAR-like dynamics persists);
BK's core ends at ~60-100 kpc (NFW-like beyond). Can current outer-halo data separate?

Models (bracketing, labeled — not a frozen-precision test):
  OURS-isolated : v from g_RAR(g_bar), a0 = 1.0801e-10 (locked), M_b fork [6,9]e10 Msun
  OURS-EFE      : boost quenched by external field e = g_ext/a0 in {0.01, 0.03}
                  via the standard 1D approximation g = g_bar * nu((g_bar + e*a0)/a0),
                  FLOORED at Newtonian-on-total (baryons + condensate halo, shared
                  M200 = 1.1e12): when the phonon boost dies, the condensate's real
                  mass still gravitates.
  BK / LCDM     : superfluid core to 80 kpc (matched inside), NFW(M200 = 1.1e12, c = 10)
                  + baryons beyond.
Anchors (robust literature enclosed masses, transcribed):
  M(<39.5 kpc) = 0.44 +/- 0.06 e12  (Watkins+19, Gaia DR2 GCs)
  M(<100 kpc)  = 0.61 +/- 0.12 e12  (Deason+21)
  M200 = 1.1 +/- 0.2 e12 at r200 ~ 219 kpc (Gaia-era consensus, Wang+20 review)
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t3_guard import ROOT, DERIVED

G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19
A0 = 1.0801e-10
nu = lambda y: 1.0/(1.0 - np.exp(-np.sqrt(y)))          # RAR interpolation
def v_kms(g, r): return np.sqrt(g*r)/1e3

r = np.geomspace(30, 320, 200)*KPC
anchors = [("M(<39.5)", 39.5, 0.44e12, 0.06e12, "Watkins+19"),
           ("M(<100)",  100.0, 0.61e12, 0.12e12, "Deason+21"),
           ("M200",     219.0, 1.10e12, 0.20e12, "Wang+20 consensus")]

M200 = 1.1e12*MSUN; c = 10.0; r200 = 219*KPC; rs = r200/c
mu = lambda x: np.log(1+x) - x/(1+x)
M_nfw = lambda rr: M200*mu(rr/rs)/mu(c)

def curves(Mb):
    gb = G*Mb*MSUN/r**2
    v_iso = v_kms(gb*nu(gb/A0), r)
    v_flo = v_kms(G*(M_nfw(r) + Mb*MSUN)/r**2, r)        # Newtonian on total (shared halo)
    out = {"iso": v_iso}
    for e in (0.01, 0.03):
        v_efe = v_kms(gb*nu(gb/A0 + e), r)
        out[f"efe{e}"] = np.maximum(v_efe, v_flo)         # mass floor when boost quenched
    out["bk"] = np.where(r < 80*KPC, v_iso, v_flo)        # BK: core matched inside, NFW+bar outside
    out["floor"] = v_flo
    return out

print("=== D24: R_SF outer-halo test ===\n--- anchor pulls (sigma) ---")
res = {}
for Mb in (6.5e10,):
    cv = curves(Mb)
    print(f"M_b = {Mb:.1e} Msun    (fork [6,9]e10 shifts RAR curves by ~ +/-4%)")
    print(f"{'anchor':10s} {'meas v':>7s}  {'ours-iso':>8s} {'ours-efe.03':>11s} {'BK/NFW':>7s}")
    for name, rk, M, dM, ref in anchors:
        vmeas = v_kms(G*M*MSUN/(rk*KPC)**2, rk*KPC); dv = 0.5*dM/M*vmeas
        i = np.argmin(np.abs(r - rk*KPC))
        row = {k: float(cv[k][i]) for k in cv}
        pulls = {k: (row[k]-vmeas)/dv for k in ("iso", "efe0.03", "bk")}
        print(f"{name:10s} {vmeas:6.0f}±{dv:.0f}  {row['iso']:7.0f}({pulls['iso']:+.1f})"
              f" {row['efe0.03']:9.0f}({pulls['efe0.03']:+.1f}) {row['bk']:6.0f}({pulls['bk']:+.1f})  [{ref}]")
        res[name] = dict(v_meas=vmeas, dv=dv, models=row, pulls=pulls)

i150 = np.argmin(np.abs(r-150*KPC)); i50 = np.argmin(np.abs(r-50*KPC)); i225 = np.argmin(np.abs(r-225*KPC))
cv = curves(6.5e10)
disc = {k: float(cv[k][i150]/cv[k][i50]) for k in ("iso", "efe0.01", "efe0.03", "bk")}
print("\n--- the discriminant: v_c(150)/v_c(50) ---")
for k, v in disc.items(): print(f"  {k:8s}: {v:.3f}")
dpc = 100*(disc["efe0.03"]-disc["bk"])
print(f"model separation (ours-EFE vs BK): {dpc:+.1f}% in the shape ratio; at 225 kpc: "
      f"{100*(cv['efe0.03'][i225]-cv['bk'][i225])/cv['bk'][i225]:+.1f}% in v_c")
print("current anchor errors: ~7-10% in v -> CANNOT separate today")
print("decisive: v_c at 150-250 kpc to <5% (Gaia-era stream modeling is approaching this)")

print("\n--- verdict ---")
print("Both models sit within ~1-2 sigma of all three anchors (ours-isolated runs high at")
print("r200, +2.3 sigma, but the EFE-quenched + mass-floor variant — the physical one at")
print("these radii — is compatible everywhere). The theories CONVERGE at large r by our")
print("own physics (boost quenches -> shared halo mass floor): the discriminant zone is")
print("100-250 kpc where the boost is active-but-marginal, and it is a 4-8% effect in v_c.")
print("STATUS: INCONCLUSIVE WITH CURRENT DATA — quantified, with the decisive dataset named.")

json.dump(dict(anchors=res, shape_ratio=disc,
               separation_pct_at_225=float(100*(cv['efe0.03'][i225]-cv['bk'][i225])/cv['bk'][i225]),
               verdict="inconclusive today; 4-8% v_c effect at 150-250 kpc; needs <5% stream-based v_c",
               notes=["ours-isolated high at r200 (+2.3 sig) but EFE+mass-floor variant compatible",
                      "convergence at large r is OUR OWN physics (quench -> condensate mass floor)",
                      "bracketing analysis, literature anchors transcribed; not a frozen-precision test"]),
          open(os.path.join(DERIVED, "b2_rsf_outer.json"), "w"), indent=1)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6.2, 3.8))
rk = r/KPC
ax.fill_between(rk, cv["efe0.03"], cv["iso"], color="crimson", alpha=.15, lw=0,
                label="ours: EFE-quenched (e=0.03) ↔ isolated")
ax.plot(rk, cv["iso"], "crimson", lw=1.2)
ax.plot(rk, cv["efe0.03"], "crimson", lw=1.2, ls=":")
ax.plot(rk, cv["bk"], "#1f77b4", lw=1.6, ls="--", label="BK: superfluid core (<80 kpc) + NFW halo")
ax.plot(rk, cv["floor"], ".5", lw=.9, ls="-.", label="shared Newtonian mass floor (M200=1.1e12)")
for name, rr, M, dM, ref in anchors:
    vmeas = np.sqrt(G*M*MSUN/(rr*KPC))/1e3; dv = 0.5*dM/M*vmeas
    ax.errorbar([rr], [vmeas], yerr=[dv], fmt="ko", ms=5, capsize=3)
    ax.annotate(ref, (rr, vmeas), xytext=(4, -14), textcoords="offset points", fontsize=6.5)
ax.axvspan(185, 310, color="crimson", alpha=.06, lw=0); ax.text(230, 236, "our $R_{SF}$", fontsize=7, color="crimson")
ax.axvspan(60, 100, color="#1f77b4", alpha=.08, lw=0); ax.text(64, 236, "BK $R_{SF}$", fontsize=7, color="#1f77b4")
ax.set(xscale="log", xlim=(30, 320), ylim=(120, 245),
       xlabel=r"$r$ [kpc]", ylabel=r"$v_c$ [km s$^{-1}$]")
ax.set_xticks([30, 50, 100, 150, 200, 300]); ax.set_xticklabels(["30", "50", "100", "150", "200", "300"])
ax.legend(loc="lower left", fontsize=7)
ax.set_title("D24: MW outer halo — the theories converge where they differ", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures/fig_rsf_outer.png"), dpi=220, bbox_inches="tight")
print("\nwrote derived/b2_rsf_outer.json, figures/fig_rsf_outer.png")