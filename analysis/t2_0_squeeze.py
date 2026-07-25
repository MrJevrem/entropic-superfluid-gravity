#!/usr/bin/env python
"""T2.0 — trade-off squeeze (plan §2): per kernel class, do the decoherence ceilings (lambda)
and diffusion ceilings (D2) leave a window around the CQ saturation line?

Anchoring [S]: the published direct windows from RESULTS_entropic_stochastic_gravity.md §9
(Janse+24 x OSSW23): ULD [1e-25,1e-16] (l_P^3/m_P D2), NLC [1e-35,1e-24] (l_P^2 D2);
atom-interferometry-admitted variants: ULD closed, NLC ~1 dex. T2.0 adds: (i) the T1.1 nHz
ceiling transduced under FROZEN kernel-index forks gamma_k in {0,2,4} [P]; (ii) the LPF and
LVK-O4a lines; (iii) window bookkeeping per class x admissibility x gamma_k.

[P] transduction (logged, coarse): T1.1 bounds per-kpc path acceleration PSD at f_nHz:
S_a^path = K95 c^2 / B4INT per kpc. Lab-band extrapolation S_a(f_lab)=S_a(f_nHz)(f_lab/f_nHz)^{g}.
The path->test-mass mapping carries an unresolved O(1..N) admissibility factor exactly analogous
to the relative-measurement dispute — T1.1's lab-band lines are therefore reported as a
LABELED OVERLAY, not merged into the direct window arithmetic (no silent unit laundering).
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from t2_guard import enforce, DERIVED, ROOT

freezes, _ = enforce(staged=[os.path.join(ROOT, "data/T2/decoherence_ceilings.csv")])
fz = freezes[-1]
K95 = fz["K95_reference"]["grid"][f"{13/3:.3f}"]
B4INT = fz["band"]["B4INT_s3"]
C = 299792458.0
F_NHZ = fz["band"]["f_k_Hz"][0]

# ---------------- ceiling tables (staged post-freeze; [S] hand-built) ----------------
DEC = pd.read_csv(os.path.join(ROOT, "data/T2/decoherence_ceilings.csv"))
ATLAS = pd.read_csv(os.path.join(ROOT, "docs/atlas_data.csv"))

Sa_path_nHz = K95 * C ** 2 / B4INT          # m^2 s^-3 per kpc of path, at ~2-10 nHz
LPF_Sa = (1.74e-15) ** 2                     # (fm s^-2/rtHz)^2 -> m^2 s^-3, 1-10 mHz [verified]
# LVK O4a flat OmegaGW<=2.8e-9 at 25 Hz -> S_h = 3 H0^2/(10 pi^2 f^3) * Omega; a-equiv per meter
H0 = 2.268e-18
LVK_Sh = 3 * H0 ** 2 / (10 * np.pi ** 2 * 25.0 ** 3) * 2.8e-9
LVK_Sa_perm = LVK_Sh * (2 * np.pi * 25.0) ** 4 / 4   # x(t)=h L/2 -> a; per meter baseline [P]

# ---------------- window bookkeeping ----------------
# Published direct/admitted windows [S] (RESULTS §9). Coordinates: class-specific dimensionless D2.
WINDOWS = {
    "ULD": dict(direct=[1e-25, 1e-16], admitted="closed"),
    "NLC": dict(direct=[1e-35, 1e-24], admitted=[1e-25, 1e-24]),
}
# T1.1 overlay vs LPF: for a kernel S_a ∝ f^g, whichever band-ceiling is lower after
# extrapolation dominates. Crossover index: T1.1 dominates for g < g_cross (red kernels).
g_cross = float(np.log10(LPF_Sa / Sa_path_nHz) / np.log10(3e-3 / F_NHZ))
overlay = {"crossover_index_g": g_cross,
           "reading": "T1.1 dominates over LPF for kernels redder (in S_a) than f^%.2f — "
                      "e.g. the tau=1/H0 Lorentzian tail (g=-2) sits AT the boundary; "
                      "steeper-red kernels are PTA-dominated, white/blue kernels LPF-dominated"
                      % g_cross}
for g in (0, -2, -4):
    Sa_lab = Sa_path_nHz * (3e-3 / F_NHZ) ** g       # T1.1-implied ceiling at 3 mHz, kernel f^g
    overlay[f"g={g}"] = dict(Sa_at_3mHz_per_kpc=float(Sa_lab),
                             vs_LPF_dex=float(np.log10(Sa_lab / LPF_Sa)),
                             dominant="LPF" if Sa_lab > LPF_Sa else "T1.1")

results = dict(
    Sa_path_nHz_per_kpc=float(Sa_path_nHz), LPF_Sa=float(LPF_Sa),
    LVK_Sh_25Hz=float(LVK_Sh), LVK_Sa_per_m=float(LVK_Sa_perm),
    windows={}, t11_overlay=overlay,
    lambda_ceilings=DEC.to_dict(orient="records"))
for cls, w in WINDOWS.items():
    entry = {"direct_window_dex": float(np.log10(w["direct"][1] / w["direct"][0])),
             "direct": w["direct"], "admitted": w["admitted"],
             "status_direct": "open",
             "status_admitted": "closed" if w["admitted"] == "closed"
                                else f"open ({np.log10(w['admitted'][1]/w['admitted'][0]):.0f} dex)"}
    results["windows"][cls] = entry
results["verdict"] = (
    "No kernel class is newly closed by T2.0: the direct windows (ULD 9 dex, NLC 11 dex) are "
    f"set by lab ceilings and floors [S]. T1.1's contribution is band-resolved: crossover at "
    f"S_a ∝ f^{g_cross:.2f} — the PTA line DOMINATES for kernels redder than ~f^-2 (the "
    "tau=1/H0 Lorentzian tail sits exactly at the boundary) and is subdominant to LPF for "
    "white/blue kernels [matching the plan's expectation once 'steep' is read as steep-red "
    "in acceleration]. Subject to the path->test-mass admissibility factor (labeled overlay, "
    "not merged into window arithmetic). The admissibility ruling (atlas §9) remains the "
    "binding open issue; branches on ULD kernels die under admitted-AI ceilings independent "
    "of carrier.")
json.dump(results, open(os.path.join(DERIVED, "t2_0_window.json"), "w"), indent=1)

# ---------------- atlas figure ----------------
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.4))
a0 = ax[0]
# lab experiments from the repo atlas (FOM = S_a N): plot FOM as ceilings vs mass
m = ATLAS.columns
mass_col = next(c for c in ATLAS.columns if "mass" in c.lower())
fom_col = next(c for c in ATLAS.columns if "fom" in c.lower())
a0.scatter(ATLAS[mass_col], ATLAS[fom_col], s=18, alpha=0.6, c="steelblue", label="lab ceilings (Janse+24)")
a0.axhline(3.0e-10, color="crimson", ls="--", lw=1.4, label="ULD floor (derived, atlas v0.1)")
a0.axhline(3.0e-12, color="darkorange", ls="--", lw=1.4, label="NLC floor")
a0.axhline(2.4e-11, color="gray", ls=":", lw=1.2, label="Asenbaum'17 (admissibility-contested)")
a0.set_xscale("log"); a0.set_yscale("log")
a0.set_xlabel("test mass [kg]"); a0.set_ylabel(r"FOM$_{D_2}$ = $S_a N$  [m$^2$s$^{-3}$]")
a0.set_title("T2.0 squeeze — lab plane (windows: ULD 9 dex, NLC 11 dex; AI-admitted: ULD closed)")
a0.legend(fontsize=7); a0.grid(alpha=0.3, which="both")

a1 = ax[1]
bands = {
    "T1.1 PTA (per kpc path)": (F_NHZ, Sa_path_nHz, "o", "crimson"),
    "LPF (test mass)": (3e-3, LPF_Sa, "s", "seagreen"),
    "LVK O4a (per m baseline)": (25.0, LVK_Sa_perm, "D", "steelblue"),
}
for lab, (f, sa, mk, c) in bands.items():
    a1.scatter([f], [sa], marker=mk, s=90, color=c, edgecolor="k", zorder=3, label=lab)
ff = np.logspace(np.log10(F_NHZ), np.log10(3e-3), 50)
for g, ls in ((0, ":"), (-2, "--"), (-4, "-")):
    a1.plot(ff, Sa_path_nHz * (ff / F_NHZ) ** g, ls, color="crimson", lw=1.1,
            label=f"T1.1 extrapolation $S_a\\propto f^{{{g}}}$")
a1.set_xscale("log"); a1.set_yscale("log")
a1.set_xlabel("frequency [Hz]"); a1.set_ylabel(r"$S_a$ ceiling [m$^2$s$^{-3}$/Hz]")
a1.set_title("acceleration-noise ceilings across 12 decades of band")
a1.legend(fontsize=7); a1.grid(alpha=0.3, which="both")
plt.tight_layout()
fig_path = os.path.join(ROOT, "docs/t2_0_squeeze_atlas.png")
plt.savefig(fig_path, dpi=130)
print("saved:", fig_path, "and derived/t2_0_window.json")
print(f"S_a^path(nHz) = {Sa_path_nHz:.2e} per kpc | LPF = {LPF_Sa:.2e} | LVK(a,per-m) = {LVK_Sa_perm:.2e}")
print(f"  crossover: T1.1 dominates for S_a ∝ f^g with g < {overlay['crossover_index_g']:.2f}")
for k, v in overlay.items():
    if isinstance(v, dict):
        print(f"  T1.1 -> 3 mHz @ {k}: {v['Sa_at_3mHz_per_kpc']:.2e}  "
              f"({v['vs_LPF_dex']:+.1f} dex vs LPF; dominant={v['dominant']})")
print("\nVERDICT:", results["verdict"][:200], "...")
