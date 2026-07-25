#!/usr/bin/env python
"""
T2 §1 + T2.1 + T2.5-phase-1 — THE FREEZE. Computes and locks ALL branch predictions BEFORE
any comparison data is staged (T2_PLAN_branch_discrimination.md).

Everything tagged [P] carries an explicit `assumptions` entry. Where a derivation forks, BOTH
forks are frozen as labeled variants (plan §11.1). Nothing in this file may be edited after
comparison data staging except via a new versioned freeze with `revision_reason`.

Contents frozen:
  T2.1  per-branch K predictions (s^2/kpc, band k=1-5 of T=16.03yr) + K95(gamma_c) reference grid
  T2.4  a0(z) locked curve (B3)
  T2.5  wide-binary prediction class + delta-v-tilde per branch (derivations logged), and the
        pipeline adjudication rules (estimator agreement, ecc priors, triple handling, thresholds)
  T2.3  early-universe amplitude scalings (B1, B5) + the thermalization dichotomy for B1
  T2.2  B4 functional prediction (achromatic-specific, slope=1 in log E vs log s_path)
"""
import os, json, subprocess, datetime
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "derived", "T2_locked_predictions.json")

# ----------------------------------------------------------------- constants (SI)
C = 299792458.0
KPC = 3.0856775814913673e19
PC = KPC / 1000
YR = 365.25 * 86400
H0 = 70.0 * 1000 / (KPC * 1000)              # 70 km/s/Mpc -> 2.268e-18 s^-1
A0 = 1.2e-10                                  # Milgrom scale, m s^-2
FYR = 1 / YR
HBAR, KB = 1.054571817e-34, 1.380649e-23

# T1.1 band (locked to the measured array bins)
T_ARR = 16.03 * YR
DF = 1 / T_ARR
FK = np.arange(1, 6) * DF                     # k=1..5: 1.98-9.88 nHz
B4INT = float(np.sum((2 * np.pi * FK) ** -4) * DF)     # s^3 — acceleration->timing band integral

def E_from_Sa(Sa_of_f):
    """Timing band power [s^2] from one-sided acceleration PSD S_a(f) [m^2 s^-4 /Hz]:
    S_r(f) = S_a/((2pi f)^4 c^2)."""
    return float(np.sum(Sa_of_f(FK) / ((2 * np.pi * FK) ** 4 * C ** 2)) * DF)

def E_from_Sr(Sr_of_f):
    return float(np.sum(Sr_of_f(FK)) * DF)

def E_powerlaw(logA, gam):
    """T1.1 v2-validated GWB-convention band power (unit test target)."""
    return float((10 ** logA) ** 2 / (12 * np.pi ** 2) * FYR ** (gam - 3)
                 * np.sum(FK ** -gam) * DF)

# ----------------------------------------------------------------- validation (plan T2.1.1)
# Push the GWB through the pipeline; ceffyl fs{hd} medians k=2..4 were -6.74,-7.14,-7.44 (T1.1 v2).
gwb_check = {f"k{k}": round(0.5 * np.log10((10 ** -14.62) ** 2 / (12 * np.pi ** 2)
             * FYR ** (13 / 3 - 3) * (k * DF) ** (-13 / 3) * DF), 3) for k in (2, 3, 4)}
VALIDATION = {"powerlaw_log10rho_k2_k4": gwb_check,
              "ceffyl_reference": {"k2": -6.74, "k3": -7.14, "k4": -7.44},
              "tolerance_dex": 0.15}

# ----------------------------------------------------------------- K95 reference grid
# Recomputed for the envelope-minimum pulsars identified in T1.1 (J1909-3744 [PPTA joint],
# J1640+2224 [NG SPNA]) across gamma_c — these set the min in v2/v3/v4; grid = their min.
def k95_grid():
    import re, glob
    res = {}
    # PPTA J1909
    PP = os.path.join(ROOT, "data/T1.1_pulsar_timing/ppta_dr3/PPTA-DR3/analysis_codes/data/all/chains")
    pars = [l.strip() for l in open(PP + "/commonNoise_pl_nocorr_freegam_DE440_pars.txt") if l.strip()]
    ch = np.load(PP + "/chain_commonNoise_pl_nocorr_freegam_DE440.npy")
    ch = ch[len(ch) // 4:]
    col = {p: i for i, p in enumerate(pars)}
    crn = ((10 ** ch[:, col["gw_pl_nocorr_freegam_log10_A"]]) ** 2 / (12 * np.pi ** 2)
           * FYR ** (ch[:, col["gw_pl_nocorr_freegam_gamma"]] - 3)
           * np.sum(FK[None, :] ** -ch[:, col["gw_pl_nocorr_freegam_gamma"]][:, None], 1) * DF)
    sources = {"J1909-3744": dict(logA=ch[:, col["J1909-3744_red_noise_log10_A"]],
                                  gam=ch[:, col["J1909-3744_red_noise_gamma"]],
                                  d=1 / 0.86, common=crn)}
    # NG J1640 (SPNA + fs-hd median common)
    ND = os.path.join(ROOT, "data/T1.1_pulsar_timing/nanograv_15yr/NANOGrav15yr_PulsarTiming_v2.1.0/narrowband/noise")
    pr = open(ND + "/J1640+2224.nb.pars.txt").read().split()
    iA = [i for i, p in enumerate(pr) if p.endswith("red_noise_log10_A")][0]
    ig = [i for i, p in enumerate(pr) if p.endswith("red_noise_gamma")][0]
    o = sorted([iA, ig])
    cc = pd.read_csv(ND + "/J1640+2224.nb.chain_1.txt", sep=r"\s+", header=None, usecols=o).values
    cm = {o[k]: cc[:, k] for k in range(2)}
    E_hd = E_powerlaw(-14.62, 13 / 3)
    sources["J1640+2224"] = dict(logA=cm[iA][len(cc) // 4:], gam=cm[ig][len(cc) // 4:],
                                 d=1 / 0.586, common=np.full(len(cc) - len(cc) // 4, E_hd))
    rng = np.random.default_rng(20260724)
    def wq(x, w, q):
        s = np.argsort(x); x, w = x[s], w[s]
        c = np.cumsum(w); c /= c[-1]
        return np.interp(q, c, x)
    for gc in (3.0, 3.5, 4.0, 13 / 3, 4.5):
        pb = (1 / (12 * np.pi ** 2)) * FYR ** (gc - 3) * np.sum(FK ** -gc) * DF
        ks = []
        for nm, s in sources.items():
            w = np.exp(-0.5 * ((s["gam"] - gc) / 0.15) ** 2)
            if w.sum() <= 0: continue
            a = s["logA"][rng.choice(len(w), 3000, p=w / w.sum())]
            E_tot = pb * 10 ** (2 * a) + s["common"][rng.integers(len(s["common"]), size=3000)]
            ks.append(float(wq(E_tot, 10 ** a, 0.95)) / s["d"])
        res[f"{gc:.3f}"] = min(ks)
    return res

K95_GRID = k95_grid()
K95 = K95_GRID[f"{13/3:.3f}"]

# ----------------------------------------------------------------- the lock (full strength)
# Program statement (Oppenheim-Russo): in the deep-MOND regime the mean entropic force EMERGES
# from the fluctuations => delta-a_rms ~ <a> ~ a0 at the field's correlation scale.
# LOCK: xi = delta a_rms / a0 = 1.  [P] — the single normalization; forks below are kernels.
XI = 1.0

# Fork COSMO (local buffeting, Lambda-coherent): tau = 1/H0 Lorentzian.
def Sa_cosmo(f):
    tau = 1 / H0
    return 4 * (XI * A0) ** 2 * tau / (1 + (2 * np.pi * f * tau) ** 2)
E_COSMO = E_from_Sa(Sa_cosmo)                 # per-pulsar-term, distance-INDEPENDENT

# Fork PATH (cell kernel, path-accumulated -> the T1.1-bounded class):
# potential cells delta_Phi = xi*a0*l_c, Shapiro fluctuation/cell dt = 2 delta_Phi l_c / c^3,
# refresh tau_c = l_c / v_eff; N = d/l_c cells; one-sided Lorentzian in f.
def K_path(l_c, v_eff):
    dt = 2 * XI * A0 * l_c ** 2 / C ** 3
    tau = l_c / v_eff
    sig2_per_kpc = (KPC / l_c) * dt ** 2
    Sr = lambda f: sig2_per_kpc * 4 * tau / (1 + (2 * np.pi * f * tau) ** 2)
    return E_from_Sr(Sr)                       # s^2 per kpc
LC_GRID = {"0.1pc": 0.1 * PC, "1pc": PC, "10pc": 10 * PC}
K_FULL = {k: K_path(l, 220e3) for k, l in LC_GRID.items()}          # matter-carrier kinematics
K_FULL_C = {k: K_path(l, C) for k, l in LC_GRID.items()}            # massless carrier (B5)

# ----------------------------------------------------------------- branch predictions
BR = {}
BR["B1"] = dict(
    carrier="all-species entanglement (induced-G)",
    K_pred={f"{lk}|p={p}": K_FULL[lk] / (107 ** p)   # evaluated AT the SM content for reference
            for lk in K_FULL for p in (0.5, 1.0)},
    N_min={f"{lk}|p={p}": float(max(1.0, (K_FULL[lk] / K95) ** (1 / p)))
           for lk in K_FULL for p in (0.5, 1.0)},
    early_universe=dict(rule="hidden dof ever thermalized above EW: DeltaNeff >= 0.027 per bosonic dof",
                        N_max_thermalized=float(0.3 / 0.027),
                        dichotomy="N_min > N_max_thermalized forces the never-thermalized loophole "
                                  "=> retreat-flag per plan §5.3"),
    wb=dict(cls="newtonian", dv=0.0,
            derivation="species count renormalizes G globally; no extra kAU-scale force"),
    assumptions=["dilution exponent p in {1/2,1} (incoherent vs coherent species sum) — both frozen",
                 "K_full from PATH kernel, xi=1, v_eff=220 km/s, l_c grid frozen",
                 "N counted in effective bosonic dof relative to visible ~107"])
BR["B2"] = dict(
    carrier="cold dark condensate (superfluid phonons)",
    K_pred="T_c-suppressed: K = K_full * (T_c/T_dS); inverse problem frozen instead",
    T_c_ceiling_over_TdS={lk: float(K95 / K_FULL[lk]) for lk in K_FULL},
    T_dS_K=float(HBAR * H0 / (2 * np.pi * KB)),
    wb=dict(cls_fork_a="mond_efe_analog", dv_fork_a=0.065,
            cls_fork_b="newtonian", dv_fork_b=0.0,
            derivation="fork a: linearize P(X)~X sqrt|X| phonon kinetic term about the dominant "
                       "galactic gradient a_ext ~ 1.6 a0 at R_sun => pair force rides external slope; "
                       "G_eff/G - 1 ~ (1/2) sqrt(a0/a_ext) ~ 0.13 => dv = sqrt(1.13)-1 ~ 0.065. "
                       "fork b: condensate coherence absent at kAU (Berezhiani-Khoury solar-system "
                       "caveat) => pure Newton."),
    assumptions=["phonon FDT amplitude ∝ k_B T_bath (linear)", "reference temperature = de Sitter T_dS",
                 "EFE-analog linearization coefficient 1/2 — order-unity uncertainty logged"])
BR["B3"] = dict(
    carrier="cosmological-horizon entropy (IR/elastic)",
    K_pred_leakage_E_per_pulsar=E_COSMO,
    leakage_margin_dex=float(np.log10(6.7e-13 / max(E_COSMO, 1e-300))),
    pta_note="distance-INDEPENDENT (local buffeting): bounded by per-pulsar quiet E95 ~ 6.7e-13 s^2, "
             "not by K95*d; f<<1/T power additionally absorbed by timing-model quadratics — "
             "null NOT over-claimable (plan T2.1.4)",
    a0_of_z=dict(rule="a0(z) = a0 * H(z)/H0, flat LCDM Om=0.3",
                 values={f"z={z}": float(np.sqrt(0.3 * (1 + z) ** 3 + 0.7))
                         for z in (0, 0.5, 1.0, 1.5, 2.0, 2.5)},
                 kill="cH(z) scaling rejected >=3sigma with inflated systematics => lock broken"),
    wb=dict(cls_fork_a="mond_no_efe_like", dv_fork_a=0.20,
            cls_fork_b="newtonian", dv_fork_b=0.0,
            derivation="fork a: Verlinde apparent-DM g_D = sqrt((cH0/6) g_B) applied to a point "
                       "mass: boost (g_D/g_B) = sqrt(a0'/g_B) — reaches ~1.7 at 30 kAU => large "
                       "rising signal, no-EFE-like. fork b: binary sits inside the Galaxy's entropy "
                       "displacement => local elastic response saturated/screened => Newtonian."),
    assumptions=["COSMO kernel Lorentzian tau=1/H0", "a0' = cH0/6 in elastic formula",
                 "screening dichotomy frozen as forks a/b — Verlinde theory does not fix it"])
BR["B4"] = dict(
    carrier="coarse-grained thermal entropy of traversed matter",
    K_pred="K(sightline) = K_full * s_path/s_ref (s_ref = sample median); functional form frozen",
    t22_prediction=dict(response="ACHROMATIC per-pulsar E_RN (v4 joint budgets)",
                        slope_logE_logS=1.0,
                        requirement="dependence must exceed the chromatic-channel dependence "
                                    "(differential statistic); loud/quiet contrast = s_max/s_min "
                                    "computed from the covariate table at comparison time",
                        exclusion_rule="measured coupling >=1 dex below lock => B4 excluded"),
    wb=dict(cls="newtonian", dv=0.0,
            derivation="kAU path entropy ~ 1e-8 of kpc sightline => force modification negligible"),
    assumptions=["linear coupling of diffusion to local entropy density (slope=1)",
                 "two-phase WIM model for s_path: T_WIM=8e3 K uniform; hot fraction from EM excess"])
BR["B5"] = dict(
    carrier="conformal/massless sector (theorem-native)",
    K_pred={f"{lk}|{sk}": K_FULL_C[lk] * sv for lk in K_FULL_C
            for sk, sv in (("Omega_r", 9.0e-5), ("Omega_r/Omega_m", 9.0e-5 / 0.31))},
    early_universe=dict(amplitude_of_z="eps(z) = eps0 * rho_r(z)/rho_tot(z) -> full strength z>3400",
                        eps0_forks=[1e-2, 1e-1],
                        observables="DeltaNeff = 7.44 * rho_fluc/rho_gamma at BBN; "
                                    "mu ~ 1.4 * d(injected)/rho_gamma over z=5e4-2e6"),
    wb=dict(cls="newtonian", dv=0.0, derivation="Omega_r-suppressed today; no kAU force"),
    assumptions=["massless carrier => v_eff = c in PATH kernel",
                 "today-suppression forks: Omega_r vs Omega_r/Omega_m — both frozen"])

# ----------------------------------------------------------------- T2.5 adjudication rules
WB_RULES = dict(
    estimators=["chae_style_deprojection", "banik_style_forward_model"],
    agreement="both must agree in sign and amplitude class for any claim",
    ecc_priors=["thermal f(e)=2e", "power-law alpha(r_p) per Hwang+2022"],
    triple_handling="RUWE<1.4 cut + nuisance f_triple in [0,0.15] marginalized",
    sample="El-Badry-Rix-Heintz eDR3; parallax/sigma>20; d<200 pc; r_p in [2,30] kAU; "
           "no tertiary flags; both components G<18",
    classes=dict(newtonian="|dv|<0.03", mond_efe="0.03<=dv<0.12", mond_no_efe="dv>=0.12"),
    threshold="class assignment requires >=3 sigma separation in the agreed estimator pair")

# ----------------------------------------------------------------- write freeze
commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"],
                        capture_output=True, text=True).stdout.strip()
dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain"],
                       capture_output=True, text=True).stdout.strip()
freeze = dict(
    provenance=dict(git_commit=commit, git_dirty=bool(dirty),
                    utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    plan="docs/T2_PLAN_branch_discrimination.md"),
    band=dict(T_yr=16.03, f_k_Hz=[float(f) for f in FK], B4INT_s3=B4INT),
    validation=VALIDATION,
    lock=dict(xi=XI, statement="delta a_rms / a0 = 1 at the field correlation scale "
              "(fluctuations generate the mean; Oppenheim-Russo mechanism)"),
    kernels=dict(cosmo="Lorentzian tau=1/H0 (local)",
                 path="Shapiro cells: l_c grid {0.1,1,10} pc, v_eff in {220 km/s, c}",
                 K_full_matter=K_FULL, K_full_massless=K_FULL_C, E_cosmo_local=E_COSMO),
    K95_reference=dict(grid=K95_GRID, source="T1.1 envelope-min pulsars J1909-3744[PPTA joint] "
                       "+ J1640+2224[NG SPNA], uniform-A 95%, recomputed per gamma_c"),
    branches=BR,
    wb_adjudication=WB_RULES,
    external_numbers_verified=dict(
        lvk_o4a="OmegaGW(25Hz) <= 2.0e-9 (a=2/3) / 2.8e-9 (flat), 95%, arXiv:2508.20721",
        neff="Neff = 3.031 +/- 0.130(stat) +/- 0.045(fg), Planck+ACT+SPT 2025 => DeltaNeff<0.3 (2sig)",
        firas="mu<9e-5, y<1.5e-5 (Fixsen+96)",
        lpf="~1.74 fm s^-2 Hz^-1/2 above 2 mHz (Armano+18)"))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(freeze, open(OUT, "w"), indent=1)
print(f"FROZEN -> {OUT}")
print(f"commit={commit[:8]} dirty={bool(dirty)}")
print(f"B4INT={B4INT:.3e} s^3 | E_cosmo={E_COSMO:.3e} s^2 | K95(13/3)={K95:.3e}")
print("K_full (matter, per l_c):", {k: f"{v:.2e}" for k, v in K_FULL.items()})
print("K_full (massless):       ", {k: f"{v:.2e}" for k, v in K_FULL_C.items()})
print("B1 N_min:", {k: f"{v:.1f}" for k, v in BR["B1"]["N_min"].items()})
print("K95 grid:", {k: f"{v:.2e}" for k, v in K95_GRID.items()})
print("validation:", VALIDATION["powerlaw_log10rho_k2_k4"], "vs", VALIDATION["ceffyl_reference"])
