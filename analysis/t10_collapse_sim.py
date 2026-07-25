#!/usr/bin/env python
"""
D29 / T-10 — the collapse simulation (toy scale; user directive 2026-07-25).

QUESTION (from D28): statics says the coherent hydrostatic state is the ground
state, with the granulation gap = 0.83x the binding energy — an O(1) photo-finish
that statics cannot call. Does violent relaxation granulate the born-coherent
field before entropic mode selection settles it?

MODEL: 2D Gross-Pitaevskii-Poisson with the quintic entropic term (hbar = m = 1):
    i dpsi/dt = [ -grad^2/2 + Phi + gamma |psi|^4 ] psi,   grad^2 Phi = A (n - nbar)
    eps(n) = gamma n^3 / 3   (the D28 convex functional; c_s^2 = 2 gamma n^2)
Born-coherent initial data (uniform phase, Gaussian blob, 1% density seed noise);
split-step spectral evolution; c_s^2 ~ n^2 turns the entropic term on LATE in the
collapse, mimicking the physical ordering (cold CDM-like infall, entropic endgame).

SCAN: gamma = 0 (validation anchor: must granulate — the known fuzzy-DM outcome)
plus gamma values targeting Theta = (5/3) c_s^2(nbar_halo)/<v^2> in [0.3 .. 2.2],
bracketing the physical point Theta_phys = 0.83. Theta is re-MEASURED per run in
the final state (mass-weighted halo density; thermal <v^2>).

PRE-REGISTERED READ-OUT (written before any run):
  per run, in the half-mass core: phase coherence C = |int psi|^2/(int sqrt(n))^2
  and granularity G = <(n - n_smooth)^2>/<n_smooth^2> (smoothing 3 lambda_dB).
  COHERENT if C > 0.7 and G < 0.3;  GRANULAR if C < 0.3 and G > 0.7;
  else PARTIAL (report radial structure: core-halo outcomes are physical).
  Verdict for the theory = state of the run(s) bracketing Theta_phys = 0.83.

Toy-scale caveats (stated up front): 2D; ~10-15 de Broglie modes across the halo
(thin speckle statistics); single isolated collapse, no mergers/expansion;
splitting error monitored via total-energy drift. A 3D + merger-history run is
the named production upgrade.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.ndimage import gaussian_filter
from t3_guard import DERIVED, ROOT

rng = np.random.default_rng(2910)
NG = 256; L = 1.0; dx = L/NG
x = (np.arange(NG) - NG/2)*dx
X, Y = np.meshgrid(x, x, indexing="ij")
k1 = 2*np.pi*np.fft.fftfreq(NG, d=dx)
KX, KY = np.meshgrid(k1, k1, indexing="ij")
K2 = KX**2 + KY**2
K2inv = np.where(K2 > 0, 1/np.maximum(K2, 1e-30), 0.0)

A_GRAV = 6.7e5          # 4 pi G in box units (sets v_vir ~ 150-190; speckle regime)
R0 = 0.18               # initial blob radius
DT = 1.0e-6; NSTEP = 16000
KIN = np.exp(-0.5j*K2*DT)

def initial_field():
    n0 = np.exp(-(X**2 + Y**2)/(2*R0**2))
    n0 *= 1.0/(n0.sum()*dx*dx)                      # total mass 1
    n0 *= (1 + 0.01*rng.standard_normal((NG, NG)))  # 1% seed noise
    return np.sqrt(np.abs(n0)).astype(np.complex128)  # uniform (zero) phase: born coherent

def poisson(n):
    rhs = np.fft.fft2(A_GRAV*(n - n.mean()))
    return np.real(np.fft.ifft2(-rhs*K2inv))

def evolve(gamma):
    psi = initial_field()
    E0 = None; drift = 0.0
    for s in range(NSTEP):
        n = np.abs(psi)**2
        V = poisson(n) + gamma*n**2
        psi *= np.exp(-0.5j*V*DT)
        psi = np.fft.ifft2(KIN*np.fft.fft2(psi))
        n = np.abs(psi)**2
        V = poisson(n) + gamma*n**2
        psi *= np.exp(-0.5j*V*DT)
        if s % 4000 == 0 or s == NSTEP-1:
            n = np.abs(psi)**2
            Ek = float(np.sum(0.5*np.abs(np.fft.ifft2(1j*KX*np.fft.fft2(psi)))**2
                              + 0.5*np.abs(np.fft.ifft2(1j*KY*np.fft.fft2(psi)))**2)*dx*dx)
            Eg = float(0.5*np.sum(poisson(n)*n)*dx*dx)
            Ee = float(np.sum(gamma*n**3/3)*dx*dx)
            Etot = Ek + Eg + Ee
            if E0 is None: E0 = Etot
            drift = abs(Etot - E0)/max(abs(E0), 1e-12)
    return psi, drift

def diagnose(psi, gamma):
    n = np.abs(psi)**2
    m_tot = n.sum()
    cx = (X*n).sum()/m_tot; cy = (Y*n).sum()/m_tot
    r = np.hypot(X - cx, Y - cy)
    # half-mass radius
    order = np.argsort(r.ravel()); csum = np.cumsum(n.ravel()[order])
    r_h = r.ravel()[order][np.searchsorted(csum, 0.5*csum[-1])]
    # thermal velocity and instantaneous de Broglie scale
    Ek = (0.5*np.abs(np.fft.ifft2(1j*KX*np.fft.fft2(psi)))**2
          + 0.5*np.abs(np.fft.ifft2(1j*KY*np.fft.fft2(psi)))**2).sum()*dx*dx
    v2 = 2*Ek/(m_tot*dx*dx)
    lam = 2*np.pi/np.sqrt(v2)
    smooth = float(np.clip(3*lam/dx, 3.0, r_h/(3*dx)))   # << r_h or granularity is meaningless
    n_s = gaussian_filter(n, smooth)
    out = {}
    for tag, sel in (("core", r < r_h), ("halo", (r >= r_h) & (r < 2.5*r_h))):
        G = float((((n[sel] - n_s[sel])**2).mean())/max((n_s[sel]**2).mean(), 1e-30))
        C = float(np.abs(psi[sel].sum())**2/max((np.sqrt(n[sel])).sum()**2, 1e-30))
        out[tag] = dict(G=round(G, 3), C=round(C, 3))
    nbar_h = float((n[r < r_h]**2).sum()/n[r < r_h].sum())      # mass-weighted density
    cs2 = 2*gamma*nbar_h**2
    theta = (5/3)*cs2/v2
    out.update(r_h=float(r_h), v_rms=float(np.sqrt(v2)), lam=float(lam),
               modes_across=float((2*r_h/lam)**2), theta=float(theta))
    return out, n

def verdict(d):
    c, g = d["core"]["C"], d["core"]["G"]
    if c > 0.7 and g < 0.3: return "COHERENT"
    if c < 0.3 and g > 0.7: return "GRANULAR"
    return "PARTIAL"

print("=== D29 / T-10: quintic-GPP collapse, Theta scan (toy scale, 2D 256^2) ===\n")
print(f"grid {NG}^2, dt = {DT:g}, steps = {NSTEP}, 4piG = {A_GRAV:g}; read-out rule pre-registered in header\n")

# pilot (gamma = 0) fixes the density scale for the Theta targets
psi0, dr0 = evolve(0.0)
d0, n0map = diagnose(psi0, 0.0)
nref = (np.abs(psi0)**2); mref = nref.sum()
cxy = ((X*nref).sum()/mref, (Y*nref).sum()/mref)
rr = np.hypot(X - cxy[0], Y - cxy[1])
order = np.argsort(rr.ravel()); csum = np.cumsum(nref.ravel()[order])
rh0 = rr.ravel()[order][np.searchsorted(csum, 0.5*csum[-1])]
nbar0 = float((nref[rr < rh0]**2).sum()/nref[rr < rh0].sum())
v20 = d0["v_rms"]**2
gate = "PASS" if d0["modes_across"] >= 8 else "FAIL — NOT in speckle regime, runs unusable"
print(f"[gamma=0 anchor] v_rms = {d0['v_rms']:.0f}, lambda = {d0['lam']:.3f}, modes across halo ~ {d0['modes_across']:.0f}, "
      f"E-drift {dr0*100:.1f}%  [speckle-regime gate: {gate}]")
print(f"                 core C = {d0['core']['C']}, G = {d0['core']['G']}  ->  {verdict(d0)} "
      f"(fuzzy-DM granulation expected: anchor {'PASS' if verdict(d0)=='GRANULAR' else 'CHECK'})\n")

results = {"anchor_gamma0": d0 | dict(verdict=verdict(d0), E_drift=dr0)}
maps = {0.0: n0map}
for th_t in (0.3, 0.6, 0.85, 1.3, 2.2):
    gamma = th_t*0.6*v20/(2*nbar0**2)        # Theta = (5/3) * 2 gamma nbar^2 / v^2
    psi, dr = evolve(gamma)
    d, nmap = diagnose(psi, gamma)
    maps[th_t] = nmap
    v = verdict(d)
    results[f"theta_target_{th_t}"] = d | dict(gamma=float(gamma), verdict=v, E_drift=dr)
    print(f"[target Th={th_t:4.2f}] measured Theta = {d['theta']:5.2f}  core: C = {d['core']['C']:.2f}, "
          f"G = {d['core']['G']:.2f}  halo: C = {d['halo']['C']:.2f}, G = {d['halo']['G']:.2f}  "
          f"drift {dr*100:.1f}%  ->  {v}")

print("\n--- figure ---")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6))
for ax, (tht, nmap) in zip(axes.ravel(), maps.items()):
    ax.imshow(np.log10(nmap.T + 1e-3), origin="lower", cmap="magma")
    key = "anchor_gamma0" if tht == 0.0 else f"theta_target_{tht}"
    d = results[key]
    th_lab = "0 (anchor)" if tht == 0.0 else f"{d['theta']:.2f}"
    ax.set_title(f"$\\Theta$={th_lab}  C={d['core']['C']:.2f} G={d['core']['G']:.2f}  {d['verdict']}",
                 fontsize=8.5)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("D29/T-10: quintic-GPP collapse — final density (log), coherence vs $\\Theta$; physical point $\\Theta$=0.83",
             fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "figures/fig_t10_collapse.png"), dpi=150)
print("wrote figures/fig_t10_collapse.png")

json.dump(results, open(os.path.join(DERIVED, "t10_collapse.json"), "w"), indent=1)
print("wrote derived/t10_collapse.json")
