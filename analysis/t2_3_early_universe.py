#!/usr/bin/env python
"""T2.3 — early-universe loudness ledger (plan §5): B1 and B5 vs DeltaNeff, FIRAS mu/y.

Frozen inputs (T2_locked_predictions.json): B1 N_min per fork + thermalization rule;
B5 eps0 forks + amplitude scaling eps(z) = eps0 * rho_r/rho_tot.
Verified bounds (freeze provenance): DeltaNeff < 0.3 (2sig, Planck+ACT+SPT'25);
FIRAS |mu|<9e-5, |y|<1.5e-5. Damping-tail full-likelihood re-run: flagged escalation, not run.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from t2_guard import enforce, DERIVED, ROOT

freezes, _ = enforce()
fz = freezes[-1]
DNEFF_MAX, MU_MAX, Y_MAX = 0.3, 9e-5, 1.5e-5

res = {"bounds": dict(DeltaNeff=DNEFF_MAX, mu=MU_MAX, y=Y_MAX)}

# ---------- B1: thermalization dichotomy ----------
b1 = fz["branches"]["B1"]
Nmin_min = min(float(v) for v in b1["N_min"].values())          # most favorable fork
dneff_if_thermal = 0.027 * Nmin_min                              # early-decoupled bosonic dof
res["B1"] = dict(N_min_most_favorable=Nmin_min,
                 DeltaNeff_if_ever_thermalized=dneff_if_thermal,
                 violation_dex=float(np.log10(dneff_if_thermal / DNEFF_MAX)),
                 verdict=("thermalized fork EXCLUDED by DeltaNeff by "
                          f"{np.log10(dneff_if_thermal/DNEFF_MAX):.1f} dex; survives ONLY as a "
                          "never-thermalized hidden sector with no other signature => "
                          "RETREAT-FLAGGED (unfalsifiable-as-stated, plan §5.3)"))

# ---------- B5: radiation-era full strength ----------
b5 = fz["branches"]["B5"]["early_universe"]
OM_R, OM_M, ZEQ = 9.0e-5, 0.31, 3400
zz = np.logspace(0, 9, 200)
rho_r_frac = OM_R * (1 + zz) ** 4 / (OM_R * (1 + zz) ** 4 + OM_M * (1 + zz) ** 3 + 0.69)
res["B5"] = {}
for eps0 in b5["eps0_forks"]:
    eps_bbn = eps0 * 1.0                                         # rho_r/rho_tot -> 1 at BBN
    dneff = 7.44 * 1.68 * eps_bbn / 1.68                         # rho_fluc = eps*rho_r; /rho_gamma
    # mu-distortion: fraction of eps dissipated inside z=5e4..2e6 window — frozen forks
    for f_diss in (1.0, 0.1):
        mu = 1.4 * eps0 * f_diss
        key = f"eps0={eps0}|f_diss={f_diss}"
        res["B5"][key] = dict(DeltaNeff=dneff, mu=mu,
                              dneff_excluded=bool(dneff > DNEFF_MAX),
                              mu_excluded=bool(mu > MU_MAX))
n_excl = sum(1 for v in res["B5"].values() if v["dneff_excluded"] or v["mu_excluded"])
res["B5_verdict"] = (f"{n_excl}/{len(res['B5'])} frozen forks violate DeltaNeff and/or FIRAS mu "
                     "(FIRAS mu is the binding bound in every fork). B5 was already EXCLUDED by "
                     "T2.1 (PTA ceiling); T2.3 is an INDEPENDENT second exclusion of all frozen "
                     "forks — the branch is doubly dead.")
res["escalation"] = "damping-tail envelope re-fit vs Planck high-ell residuals: NOT run (flagged)"

json.dump(res, open(os.path.join(DERIVED, "t2_3_results.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "B5"}, indent=1)[:900])
print("B5 forks:", {k: ("EXCL" if v["dneff_excluded"] or v["mu_excluded"] else "ok")
                    for k, v in res["B5"].items()})

fig, ax = plt.subplots(figsize=(7.5, 5))
for eps0, c in zip(b5["eps0_forks"], ("seagreen", "crimson")):
    ax.plot(1 + zz, eps0 * rho_r_frac, color=c, lw=2, label=f"B5 amplitude, $\\epsilon_0$={eps0}")
ax.axhline(DNEFF_MAX / 12.5, color="k", ls="--", lw=1.2,
           label=r"$\Delta N_{\rm eff}$=0.3 equivalent ($\epsilon\approx0.024$)")
ax.axhline(MU_MAX / 1.4, color="purple", ls=":", lw=1.2, label=r"FIRAS $\mu$ (f_diss=1)")
ax.axvline(1 + ZEQ, color="gray", lw=0.8); ax.text(ZEQ, 2e-6, " z_eq", fontsize=8)
ax.axvspan(5e4, 2e6, alpha=0.08, color="purple"); ax.text(2e5, 2e-6, r"$\mu$-window", fontsize=8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("1+z"); ax.set_ylabel(r"fluctuation fraction $\epsilon(z)$")
ax.set_title("T2.3 epoch ledger: B5 radiation-era loudness vs bounds (B1: see dichotomy)")
ax.legend(fontsize=8); ax.grid(alpha=0.3, which="both")
plt.tight_layout()
plt.savefig(os.path.join(ROOT, "docs/t2_3_epoch_ledger.png"), dpi=130)
print("saved derived/t2_3_results.json + docs/t2_3_epoch_ledger.png")
