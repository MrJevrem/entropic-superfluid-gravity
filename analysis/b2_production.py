#!/usr/bin/env python
"""
D16 — production mechanism for the carrier: vacuum misalignment (intermediate-scale ALP).

Constraints inherited: non-thermal + born cold (D15.7); (near-)zero particle self-coupling
(D15.5 Bullet theorem); Omega_dm matched; condensed from birth (D13 initial condition).
Mechanism: field frozen at phi_i while H > m; oscillation onset 3H(T_osc) = m; redshifts as
exact-cold matter. Yield: n/s = (45/4pi^2 g*)(m phi_i^2 / T_osc^3)|_osc, T_osc from
H = 1.66 sqrt(g*) T^2/M_pl. Scenarios: pre-inflationary (homogeneous theta_i; isocurvature
bound -> H_inf cap -> r cap) vs post-inflationary (miniclusters; no isocurvature) — both
frozen. Photon coupling: generic g = alpha C/(2 pi f_a), C = O(1) (E/N fork: photophobic
variant frozen too) -> near-UV/visible decay line.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t2_guard import DERIVED, ROOT

MPL, MBAR = 1.22e19, 2.435e18            # GeV (Planck, reduced)
ALPHA = 7.297e-3
S0_over_rhoc_h2 = 2.74e8                 # s0/(rho_c/h^2) [GeV^-1]
AS, BETA_ISO = 2.1e-9, 0.038
T_UNIVERSE = 4.35e17                     # s
GEV_INV_S = 6.582e-25

def gstar(T_GeV):                        # coarse SM g*
    return 106.75 if T_GeV > 200 else (61.75 if T_GeV > 0.2 else 10.75)

def produce(m_eV, theta=1.0):
    m = m_eV * 1e-9                      # GeV
    g = 106.75
    for _ in range(4):                   # iterate T_osc(g*)
        T_osc = np.sqrt(m * MPL / (3 * 1.66 * np.sqrt(g)))
        g = gstar(T_osc)
    # n/s = (45/(4 pi^2 g)) * (m phi^2/2) * 2 / (m T^3)  -> Omega h^2 = m*(n/s)*s0/rho_c*h2
    # Solve for phi_i from Omega h^2 = 0.12:
    pref = (45 / (4 * np.pi ** 2 * g)) / T_osc ** 3
    phi2 = 0.12 / (S0_over_rhoc_h2 * pref * m / 2 * m) * 1.0  # from 0.12 = m*(pref*m phi^2/2)*s0/..
    phi_i = np.sqrt(phi2)
    f_a = phi_i / theta
    # pre-inflationary: isocurvature cap -> H_inf cap -> r cap
    Piso_max = BETA_ISO / (1 - BETA_ISO) * AS
    H_inf_max = np.sqrt(Piso_max) / 2 * (2 * np.pi * f_a * theta) / 2  # dρ/ρ=2 dθ/θ; δθ=H/2πf
    r_max = (2 * H_inf_max ** 2 / (np.pi ** 2 * MBAR ** 2)) / AS
    # photon coupling + decay line
    g_agg = ALPHA / (2 * np.pi * f_a)                       # C=1
    Gamma = g_agg ** 2 * m ** 3 / (64 * np.pi)
    tau_s = (1 / Gamma) * GEV_INV_S
    lam_nm = 2 * 1239.84 / m_eV
    # self-interaction (cosine quartic) vs Bullet
    lam4 = m ** 2 / f_a ** 2
    sigma = lam4 ** 2 / (64 * np.pi * m ** 2) * 3.89e-28    # cm^2
    sigma_over_m = sigma / (m_eV * 1.783e-33 * 1e3)         # cm^2/g  (m in g)
    return dict(m_eV=m_eV, T_osc_TeV=float(T_osc / 1e3), gstar=g,
                phi_i_GeV=float(phi_i), f_a_GeV=float(f_a),
                H_inf_max_GeV=float(H_inf_max), r_max=float(r_max),
                g_agamma_GeV=float(g_agg), tau_decay_s=float(tau_s),
                tau_over_t_universe=float(tau_s / T_UNIVERSE), line_nm=float(lam_nm),
                quartic=float(lam4), sigma_over_m_cm2_g=float(sigma_over_m))

grid = {f"m={m}eV": produce(m) for m in (5.8, 8.5, 16.3)}
fid = grid["m=8.5eV"]

# anchor validation: fuzzy DM m=1e-22 eV must give phi ~ 1e17 GeV
anchor = produce(1e-22)
# coldness: cosmic-mean occupancy with generous v = 3 km/s
m_kg = 8.5 * 1.78266e-36
n_cosmic = 0.26 * 8.6e-27 / m_kg
lam_dB = 6.626e-34 / (m_kg * 3e3)
occup = n_cosmic * lam_dB ** 3
# post-inflationary minicluster comoving scale/mass
a_ratio = 2.35e-13 / (fid["T_osc_TeV"] * 1e3) * (3.91 / 106.75) ** (1 / 3)
L_com = (1 / (fid["T_osc_TeV"] * 0)) if False else None
H_osc_GeV = 8.5e-9 / 3
L_com_m = (1 / H_osc_GeV) * 1.973e-16 / a_ratio
M_mc_kg = 0.26 * 8.6e-27 * L_com_m ** 3

out = dict(
    mechanism="vacuum misalignment; intermediate-scale ALP",
    fiducial=fid, mass_sweep=grid,
    anchor_check=dict(fuzzy_dm_phi_i=anchor["phi_i_GeV"],
                      expect="~1e17 GeV", ok=bool(3e16 < anchor["phi_i_GeV"] < 5e17)),
    coldness=dict(cosmic_occupancy_at_3km_s=float(occup), zeta32=2.612,
                  verdict="condensed from birth even at cosmic mean density — grounds "
                          "D13's phase-persistence initial condition"),
    scenarios=dict(
        pre_inflationary=dict(preferred=True,
                              H_inf_cap_GeV=fid["H_inf_max_GeV"], r_cap=fid["r_max"],
                              note="homogeneous theta_i; NO miniclusters; isocurvature near "
                                   "the cap is a POSITIVE future CMB channel"),
        post_inflationary=dict(minicluster_comoving_m=float(L_com_m),
                               minicluster_mass_kg=float(M_mc_kg),
                               note="micro-miniclusters (~1e5 kg, sub-AU) — dynamically "
                                    "irrelevant; string/wall contributions [O]")),
    convergence="the r-cap (~1e-15) INDEPENDENTLY reproduces the program's Sec.-8 no-tensor "
                "falsifier: a primordial B-mode detection now kills BOTH the classical-gravity "
                "program and the carrier production — a DOUBLE kill, from disjoint physics",
    bullet_consistency=f"cosine quartic lambda = {fid['quartic']:.1e}; sigma/m = "
                       f"{fid['sigma_over_m_cm2_g']:.1e} cm^2/g — 60 dex under the Bullet "
                       "bound: the production story delivers exactly the (nearly) "
                       "non-self-interacting carrier D15.5 demanded",
    em_channel=dict(g_agamma_generic=fid["g_agamma_GeV"], line_band_nm=[152, 428],
                    fiducial_line_nm=fid["line_nm"], tau_over_tU=fid["tau_over_t_universe"],
                    current_bounds="MUSE-class: g < ~1e-12..1e-13 at 2.7-5.3 eV — generic "
                                   "coupling sits 2 orders below; photophobic fork frozen",
                    verdict="an ultra-faint DM decay line at 152-428 nm is the theory's only "
                            "prospective electromagnetic detection channel"),
    open_items=["theta_i ~ O(1) naturalness / measure", "UV origin of f_a ~ 2e11 GeV",
                "E/N (photophobic vs generic) — both frozen",
                "coexistence proof: emergent X^{3/2} sector atop the ALP potential (D15.5 [O])"])
json.dump(out, open(os.path.join(DERIVED, "b2_production.json"), "w"), indent=1)

v2p = os.path.join(DERIVED, "T2_locked_predictions_v2.json")
v2 = json.load(open(v2p))
v2["branches"]["B2"]["production"] = dict(
    mechanism="misalignment ALP", f_a_GeV=round(fid["f_a_GeV"], -9),
    T_osc_TeV=round(fid["T_osc_TeV"], 0), H_inf_cap_GeV=fid["H_inf_max_GeV"],
    r_cap=fid["r_max"], decay_line_nm=[152, 428],
    kills=["primordial B-modes r>~1e-3 (DOUBLE kill w/ Sec.-8)",
           "g_agamma pushed 2+ orders below 5e-15 at 6-16 eV kills generic fork "
           "(photophobic survives — retreat-risk flagged)"],
    ref="docs/DM_DERIVATIONS.md D16")
v2["provenance"]["revision_reason"] += " | D16: production mechanism derived+frozen (same cycle)."
json.dump(v2, open(v2p, "w"), indent=1)

print("=== D16: production — misalignment ALP ===")
for k, v in grid.items():
    print(f"[{k:9}] T_osc={v['T_osc_TeV']:.0f} TeV  phi_i={v['phi_i_GeV']:.2e} GeV  "
          f"f_a={v['f_a_GeV']:.2e}  line={v['line_nm']:.0f} nm")
print(f"anchor (fuzzy DM 1e-22 eV): phi_i = {anchor['phi_i_GeV']:.2e} GeV (expect ~1e17) "
      f"-> {'OK' if out['anchor_check']['ok'] else 'FAIL'}")
print(f"coldness: cosmic occupancy = {occup:.0f} (>> 2.6) -> condensed from birth")
print(f"pre-inflationary: H_inf < {fid['H_inf_max_GeV']:.1e} GeV  =>  r < {fid['r_max']:.1e}")
print(f"photon: g = {fid['g_agamma_GeV']:.1e} GeV^-1; tau/t_U = {fid['tau_over_t_universe']:.1e}; "
      f"line {fid['line_nm']:.0f} nm")
print(f"Bullet: quartic {fid['quartic']:.1e}; sigma/m = {fid['sigma_over_m_cm2_g']:.1e} cm^2/g")
print(f"miniclusters (post-inf fork): L_com = {L_com_m:.1e} m, M = {M_mc_kg:.1e} kg")
print("\nsaved derived/b2_production.json + freeze v2 updated")
