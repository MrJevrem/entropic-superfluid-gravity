#!/usr/bin/env python
"""
D31 — 3D dissipative vortex-tangle decay in the quintic medium (T-3 keystone;
user directive 2026-07-25).

CORRECTION TO D30'S REQUIREMENT: the ordered end state is NOT a vortex-free
30-kpc phase (the halo rotates) — it is the D23 lattice (spacing ~0.02 AU).
Phase ordering is therefore VORTEX-SPACING COARSENING delta(t), with the quantum
of circulation kappa = h/m = 44 m^2/s (D23) as the natural diffusivity, and the
Vinen-class law delta^2 ~ nu' t with nu' = O(0.1-1) kappa [K: T=0 helium
(Walmsley-Golov) and GPE turbulence simulations].

REQUIREMENT LADDER (each level: what the phenomenology needs, when it is reached):
  R1  granule DENSITY contrast decay (D28 selection operative)      ~ lambda/c_s
  R2  tangle velocity < 1% of c_s at kpc scales (P(X) hydrodynamics
      clean; residual tangle a <1e-4 energy perturbation)           delta > kappa/(2 pi 0.01 c_s)
  R3  approach to the D23 equilibrium lattice (full order)          delta -> delta_eq

Part B: 3D damped quintic GPE (64^3), speckle initial data at Theta ~ 0.85,
Lambda = 0 (T=0) and 0.05 (mutual-friction proxy): measures (i) contrast decay,
(ii) line-density decay -> the quintic medium's nu'/kappa — the one number that
anchors the halo-scale rates. PASS criterion (pre-stated): nu'/kappa >= 0.05.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.ndimage import gaussian_filter
from t3_guard import DERIVED

EV = 1.602176634e-19; C = 2.99792458e8; HBAR = 1.054571817e-34; H = 6.62607015e-34
PC = 3.0857e16; GYR = 3.156e16; AU = 1.496e11
m = 8.45*EV/C**2; vf = 156e3; cs = vf/np.sqrt(2)
kappa = H/m
lam = H/(m*cs); xi = HBAR/(m*cs)
T_H = 13.8*GYR
d_eq = 2.9e9                      # D23 lattice spacing (0.015-0.026 AU, mid)

print("=== D31: vortex-tangle decay — requirement ladder + 3D quintic measurement ===\n")
print("--- (A) the analytic ladder (kappa = h/m = %.1f m^2/s [D23]) ---" % kappa)
t_R1 = lam/cs
d_R2 = kappa/(2*np.pi*0.01*cs)
print(f"R1 density-contrast decay: t ~ lambda_dB/c_s = {t_R1*1e9:.1f} ns  ->  INSTANT")
for nu_frac in (0.1, 1.0):
    nu = nu_frac*kappa
    t_R2 = d_R2**2/nu
    t_R3 = d_eq**2/nu
    d_tH = np.sqrt(nu*T_H)
    print(f"nu' = {nu_frac:.1f} kappa: R2 (delta > {d_R2*1e3:.1f} mm): t = {t_R2*1e6:.0f} us;  "
          f"R3 (delta -> {d_eq/AU:.3f} AU): t = {t_R3/GYR:.0f} Gyr;  delta(t_H)/delta_eq = {d_tH/d_eq:.2f}")
print("READING: R1 and R2 — everything the kpc P(X) phenomenology needs — complete in")
print("nanoseconds to microseconds EVEN ON THE DIFFUSIVE LAW. R3 (the full lattice) takes")
print("6-60 Gyr: halos today sit between ~1x and ~4x the equilibrium line density — a")
print("standing micro-prediction (slightly enhanced tangle in young/recently-merged halos),")
print("not a failure. D30's 16-25-dex fork was an artifact of demanding a vortex-free")
print("30-kpc phase; the rotating end state dissolves it — CONDITIONAL on nu'/kappa >= 0.05")
print("holding in the QUINTIC medium, which Part B now measures.\n")

print("--- (B) 3D damped quintic GPE, 64^3: measuring nu'/kappa ---")
NG = 64; dx = 1.0
k1 = 2*np.pi*np.fft.fftfreq(NG, d=dx)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
K2 = KX**2 + KY**2 + KZ**2
rng = np.random.default_rng(31)
k0 = 2*np.pi/6.0                              # speckle scale lambda = 6 dx
DT = 0.06; NSTEP = 6000; SNAP = 400
KAPPA_SIM = 2*np.pi                           # h/m with hbar=m=1

def init_speckle():
    amp = np.exp(-(np.sqrt(K2) - k0)**2/(2*(0.25*k0)**2))
    psi = np.fft.ifftn(amp*np.exp(2j*np.pi*rng.random((NG, NG, NG))))
    return (psi/np.sqrt(np.mean(np.abs(psi)**2))).astype(np.complex128)

def line_density(psi):
    th = np.angle(psi); w = lambda a: (a + np.pi) % (2*np.pi) - np.pi
    tot = 0
    for ax1, ax2 in ((0, 1), (0, 2), (1, 2)):
        t1 = np.roll(th, -1, ax1); t2 = np.roll(t1, -1, ax2); t3 = np.roll(th, -1, ax2)
        circ = w(t1 - th) + w(t2 - t1) + w(t3 - t2) + w(th - t3)
        tot += np.count_nonzero(np.abs(circ) > np.pi)
    return tot*dx/NG**3                        # line length per volume

def contrast(psi):
    n = np.abs(psi)**2
    n_s = gaussian_filter(n, 2*6.0/2.355, mode="wrap")
    return float(np.mean((n - n_s)**2)/np.mean(n_s**2))

results = {}
for LAM in (0.0, 0.05):
    psi = init_speckle()
    n = np.abs(psi)**2
    Ek = float(np.sum(0.5*K2*np.abs(np.fft.fftn(psi)/NG**1.5)**2))/np.mean(n)/NG**0  # per-particle scale
    v2 = 2*Ek/NG**3*NG**3; v2 = float(np.mean(0.5*K2*np.abs(np.fft.fftn(psi))**2/NG**3/np.sum(np.abs(psi)**2)*2*NG**3))
    # robust v2: 2*E_kin/M
    Ekin = float(np.sum(0.5*K2*np.abs(np.fft.fftn(psi))**2)/NG**3)
    M = float(np.sum(np.abs(psi)**2))
    v2 = 2*Ekin/M
    gamma = 0.85*0.6*v2/(2*np.mean(n)**2)      # Theta = 0.85 target
    N0 = M
    KINF = np.exp((-1j - LAM)*0.5*K2*DT)
    series = []
    for s in range(NSTEP):
        n = np.abs(psi)**2
        mu_c = gamma*np.mean(n**2)/np.mean(n)  # keeps damped run near fixed density
        psi *= np.exp((-1j - LAM)*(gamma*n**2 - mu_c)*0.5*DT)
        psi = np.fft.ifftn(KINF*np.fft.fftn(psi))
        n = np.abs(psi)**2
        psi *= np.exp((-1j - LAM)*(gamma*n**2 - mu_c)*0.5*DT)
        if LAM > 0:
            psi *= np.sqrt(N0/np.sum(np.abs(psi)**2))
        if s % SNAP == 0 or s == NSTEP-1:
            L = line_density(psi)
            series.append(dict(t=round(s*DT, 1), L=float(L),
                               delta=float(1/np.sqrt(max(L, 1e-12))), G=round(contrast(psi), 3)))
    ts = np.array([p["t"] for p in series[len(series)//3:]])
    d2 = np.array([p["delta"]**2 for p in series[len(series)//3:]])
    slope = float(np.polyfit(ts, d2, 1)[0])
    nu_over_kappa = slope/KAPPA_SIM
    results[f"Lambda_{LAM}"] = dict(series=series, nu_over_kappa=float(nu_over_kappa),
                                    G_first=series[0]["G"], G_last=series[-1]["G"])
    print(f"[Lambda = {LAM:.2f}] contrast G: {series[0]['G']:.2f} -> {series[-1]['G']:.2f};  "
          f"delta: {series[0]['delta']:.1f} -> {series[-1]['delta']:.1f} dx;  "
          f"fit d(delta^2)/dt = {slope:.3f}  ->  nu'/kappa = {nu_over_kappa:.3f}")

print("\n--- (C) verdict ---")
nu_min = min(r["nu_over_kappa"] for r in results.values())
ok = nu_min >= 0.05
print(f"measured nu'/kappa in the quintic medium: "
      + ", ".join(f"{k.split('_')[1]}: {r['nu_over_kappa']:.3f}" for k, r in results.items()))
print(f"PASS criterion nu'/kappa >= 0.05: {'PASS' if ok else 'FAIL'}")
if ok:
    print("=> the halo-scale rates of Part A hold: R1/R2 (everything the force needs) complete")
    print("   in ns-us; R3 approaches the D23 lattice within O(t_H), leaving a <~4x residual")
    print("   line-density excess as a micro-prediction. THE D30 FORK IS DISSOLVED: even")
    print("   diffusive-class ordering suffices once the requirement is stated against the")
    print("   correct (rotating-lattice) end state. T-3's pass/fail criterion: PASSED.")
else:
    print("=> the quintic medium orders anomalously slowly; the D30 fork stands. RECORD AS ADVERSE.")

json.dump(dict(kappa_m2s=float(kappa),
               ladder=dict(R1_ns=float(t_R1*1e9), R2_delta_mm=float(d_R2*1e3),
                           R3_Gyr_band=[float(d_eq**2/kappa/GYR), float(d_eq**2/(0.1*kappa)/GYR)],
                           delta_tH_over_deq_band=[float(np.sqrt(0.1*kappa*T_H)/d_eq), float(np.sqrt(kappa*T_H)/d_eq)]),
               sim={k: dict(nu_over_kappa=r["nu_over_kappa"], G_first=r["G_first"], G_last=r["G_last"])
                    for k, r in results.items()},
               pass_criterion="nu'/kappa >= 0.05", passed=bool(ok),
               verdict=("D30 fork dissolved: requirement restated against the rotating-lattice end state; "
                        "R1/R2 in ns-us on the diffusive law itself; R3 in 6-60 Gyr with a <~4x residual-tangle "
                        "micro-prediction" if ok else "adverse: anomalously slow quintic ordering")),
          open(os.path.join(DERIVED, "b2_tangle_decay.json"), "w"), indent=1)
print("\nwrote derived/b2_tangle_decay.json")
