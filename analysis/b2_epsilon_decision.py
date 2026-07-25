#!/usr/bin/env python
"""
D8 — deciding the a0 coefficient: 1/6 vs 1/2pi (B2-conditional, user directive).

Three levers, all framework-internal:
 L1 (Lambda purity): the lock must be a pure-Lambda statement (D2 corollary; enforced
     empirically by our T2.4 const-preference at 3.2-6.3 sigma). Candidates must be written
     Lambda-natively; any form needing today's TOTAL H0 implies a0 ~ H(z) evolution => dead.
 L2 (magnitude, parameter-free): compare a0_pred(Lambda_Planck) against our SPARC fit.
 L3 (pedigree): the theorem supplies exactly two constants — the modular 2pi (Eq. 4) and the
     coupling alpha = 8pi (Eq. 28). A fluctuation-amplitude (second-moment) origin of the
     force scales as alpha^(-1/2) => kappa_dS/sqrt(8pi); a modular-frequency origin gives
     kappa_dS/2pi.
Candidates (Lambda-native):
   MOD:  a0 = c^2 sqrt(Lambda/3) / (2pi)          [modular frequency of the dS horizon]
   CPL:  a0 = c^2 sqrt(Lambda/(24 pi))            [= kappa_dS/sqrt(8pi); "1/6" at OmL~0.7]
   (bare 1/2pi * cH0 with TOTAL H0: excluded by L1 — needs a0 ~ H(z), killed by T2.4.)
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from t2_guard import DERIVED, ROOT

C = 299792458.0
MPC = 3.0856775814913673e22
A0, A0_STAT = 1.128e-10, 0.019e-10          # our SPARC fit (T2.4)
A0_SYS = 0.20 * A0                           # M/L normalization systematic (RAR literature ~20%)

def H_SI(h): return h * 1e3 / MPC

# ---------------- L2: parameter-free pure-Lambda comparison (Planck 2018) ----------------
H0_P, OML_P = 67.36, 0.6847
LAM = 3 * OML_P * H_SI(H0_P) ** 2 / C ** 2
a0_MOD = C ** 2 * np.sqrt(LAM / 3) / (2 * np.pi)
a0_CPL = C ** 2 * np.sqrt(LAM / (24 * np.pi))
dev = lambda p: (A0 - p) / p
sig = lambda p: (A0 - p) / np.sqrt(A0_STAT ** 2 + A0_SYS ** 2)
res = dict(
    inputs=dict(a0_sparc=A0, stat=A0_STAT, sys=A0_SYS, H0_planck=H0_P, OmL=OML_P, Lambda=LAM),
    candidates=dict(
        MOD=dict(form="c^2 sqrt(Lambda/3)/(2pi)", a0_pred=a0_MOD,
                 dev_pct=100 * dev(a0_MOD), sigma_with_sys=float(sig(a0_MOD))),
        CPL=dict(form="c^2 sqrt(Lambda/(24pi)) = kappa_dS/sqrt(8pi)", a0_pred=a0_CPL,
                 dev_pct=100 * dev(a0_CPL), sigma_with_sys=float(sig(a0_CPL))),
        ratio_CPL_over_MOD=float(a0_CPL / a0_MOD)),      # = sqrt(pi/2) exactly
)

# ---------------- identity + implied parameters ----------------
res["identities"] = dict(
    sqrt_OmL_over_8pi_at_0p7=float(np.sqrt(0.7 / (8 * np.pi))),
    one_sixth=1 / 6,
    reading="sqrt(OmL/8pi) = 0.16689 at OmL=0.7 vs 1/6 = 0.16667 (0.13%): Verlinde's 1/6 is "
            "the OmL~0.7 shadow of sqrt(1/8pi); the two differ ONLY through OmL-scaling.")
eps_cpl = np.sqrt(0.7 / (8 * np.pi))
H0_implied = A0 / (C * eps_cpl) * MPC / 1e3
res["implied"] = dict(
    H0_if_CPL=float(H0_implied),
    Lambda_from_a0=float(24 * np.pi * A0 ** 2 / C ** 4),
    Lambda_planck=LAM,
    lambda_agreement_pct=float(100 * (24 * np.pi * A0 ** 2 / C ** 4 / LAM - 1)))

# ---------------- verdict ----------------
res["verdict"] = (
    f"DECIDED FOR THE 1/6-CLASS, with derived exact form a0 = c^2 sqrt(Lambda/24pi) = "
    f"kappa_dS/sqrt(8pi). Grounds: (L1) the bare-1/2pi reading fits only at H0~73 with a lock "
    f"to TOTAL H(z=0), which implies a0~H(z) evolution — killed by our T2.4 (const preferred "
    f"3.2-6.3 sigma). (L2) parameter-free vs Planck Lambda: CPL deviates "
    f"{100*dev(a0_CPL):+.1f} pct ({sig(a0_CPL):+.2f} sigma incl. M/L systematics); MOD deviates "
    f"{100*dev(a0_MOD):+.1f} pct ({sig(a0_MOD):+.2f} sigma) — the Lambda-native modular route is "
    f"disfavored, the coupling route agrees. (L3) the sqrt(8pi) is the square root of the "
    f"paper's Eq.-(28) coupling alpha=8pi — the natural amplitude scaling for a force that is a "
    f"SECOND-MOMENT (mean-from-noise) effect, per the Oppenheim-Russo mechanism; the modular "
    f"2pi belongs to frequency conversion, not force amplitude. NOT yet a >3-sigma magnitude "
    f"exclusion of MOD on a0 alone (M/L systematic ~20 pct); named escalation: M/L-marginalized "
    f"RAR refit. Corollaries: implied H0 = {H0_implied:.1f} km/s/Mpc (sides with ~70, against "
    f"73 — the 1/2pi@H0=73 reading is doubly rejected); a0 becomes a Lambda-meter: Lambda(a0) "
    f"agrees with Planck to {res['implied']['lambda_agreement_pct']:+.0f} pct. Discriminator vs "
    f"exact-elastic-1/6: our coefficient scales as sqrt(OmL); Verlinde's is OmL-independent — "
    f"separable in principle by joint (a0, OmL, H0) fits.")
res["tier2_kill_test"] = (
    "Compute the B-coefficient of the Sec.-4 template on the acoustic horizon algebra "
    "(the 'constants bookkeeping against the 8pi normalization' the program itemizes): the "
    "identification DIES unless the Lambda-channel variance bookkeeping returns 1/8pi.")

json.dump(res, open(os.path.join(DERIVED, "b2_epsilon.json"), "w"), indent=1)

# ---------------- figure ----------------
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5))
h = np.linspace(65, 75, 200)
eps_meas = A0 / (C * H_SI(h))
a0 = ax[0]
a0.plot(h, eps_meas, "k-", lw=2, label=r"measured $\epsilon = a_0^{\rm SPARC}/cH_0$")
a0.fill_between(h, (A0 - A0_STAT) / (C * H_SI(h)), (A0 + A0_STAT) / (C * H_SI(h)),
                color="k", alpha=0.25, label="stat")
a0.fill_between(h, (A0 - A0_SYS) / (C * H_SI(h)), (A0 + A0_SYS) / (C * H_SI(h)),
                color="k", alpha=0.08, label="M/L sys")
for val, lab, c in ((1 / 6, r"$1/6$", "seagreen"),
                    (np.sqrt(0.7 / (8 * np.pi)), r"$\sqrt{\Omega_\Lambda/8\pi}$", "crimson"),
                    (1 / (2 * np.pi), r"$1/2\pi$ (needs total $H_0$: L1-excluded)", "steelblue"),
                    (np.sqrt(0.7) / (2 * np.pi), r"$\sqrt{\Omega_\Lambda}/2\pi$ (modular)", "purple")):
    a0.axhline(val, color=c, ls="--", lw=1.4, label=f"{lab} = {val:.4f}")
a0.axvline(H0_implied, color="crimson", lw=0.8, alpha=0.6)
a0.text(H0_implied + 0.1, 0.145, f"implied $H_0$={H0_implied:.1f}", fontsize=8, color="crimson")
a0.set_xlabel(r"$H_0$ [km/s/Mpc]"); a0.set_ylabel(r"$\epsilon$")
a0.set_title("D8: coefficient candidates vs measurement across the $H_0$ tension")
a0.legend(fontsize=7.5); a0.grid(alpha=0.3)

a1 = ax[1]
bars = {"CPL:  $c^2\\sqrt{\\Lambda/24\\pi}$": a0_CPL,
        "MOD:  $c^2\\sqrt{\\Lambda/3}/2\\pi$": a0_MOD}
xx = np.arange(len(bars))
a1.bar(xx, list(bars.values()), 0.5, color=["crimson", "purple"], alpha=0.8)
a1.axhline(A0, color="k", lw=2, label=f"SPARC fit {A0:.3e}")
a1.axhspan(A0 - A0_SYS, A0 + A0_SYS, color="k", alpha=0.08, label="M/L sys band")
a1.set_xticks(xx); a1.set_xticklabels(list(bars.keys()), fontsize=9)
for i, (k, v) in enumerate(bars.items()):
    a1.text(i, v * 1.01, f"{v:.3e}\n({100*(A0-v)/v:+.1f}%)", ha="center", fontsize=8)
a1.set_ylabel(r"$a_0$ [m s$^{-2}$]")
a1.set_title(r"parameter-free (Planck $\Lambda$): the $\sqrt{\pi/2}$ gap decides")
a1.legend(fontsize=8); a1.grid(alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig(os.path.join(ROOT, "docs/b2_epsilon_decision.png"), dpi=130)
print(f"MOD: a0 = {a0_MOD:.4e}  ({100*dev(a0_MOD):+.1f}%, {sig(a0_MOD):+.2f} sig)")
print(f"CPL: a0 = {a0_CPL:.4e}  ({100*dev(a0_CPL):+.1f}%, {sig(a0_CPL):+.2f} sig)")
print(f"ratio CPL/MOD = {a0_CPL/a0_MOD:.4f} (sqrt(pi/2)={np.sqrt(np.pi/2):.4f})")
print(f"sqrt(OmL/8pi)@0.7 = {eps_cpl:.5f} vs 1/6 = {1/6:.5f}")
print(f"implied H0 = {H0_implied:.2f} | Lambda(a0)/Lambda_Planck - 1 = "
      f"{res['implied']['lambda_agreement_pct']:+.0f}%")
print("saved derived/b2_epsilon.json + docs/b2_epsilon_decision.png")
