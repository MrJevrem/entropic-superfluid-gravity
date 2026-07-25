#!/usr/bin/env python
"""
D23 — vortices in the rotating entropic superfluid halo (BK-comparison open item 2).
 (1) Feynman array: circulation quantum kappa = h/m; areal density n_v = 2 Omega/kappa
     for halo spin (lambda_spin fork 0.02-0.06); spacing, total number, energy budget,
     density-deficit observability bound. BK-contrast at their m ~ 0.6 eV.
 (2) NEW (extends [N11]): the quintic vortex profile — shooting solution of
     u'' + u'/r - u/r^2 = 2(u^5 - u)   (equation units; r_xi = r*sqrt(2))
     vs the cubic vortex u'' + u'/r - u/r^2 = 2(u^3 - u); core widths at matched xi.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from t3_guard import ROOT, DERIVED

H = 6.62607015e-34; EV = 1.602176634e-19; C = 2.99792458e8
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19
m = 8.45 * EV / C**2; m_bk = 0.6 * EV / C**2
vf = 220e3; cs = vf/np.sqrt(2); xi = (H/(2*np.pi))/(m*cs)
RSF = 200*KPC

print("=== D23: vortices in the rotating halo ===\n--- (1) the Feynman array ---")
kappa = H/m; kappa_bk = H/m_bk
print(f"circulation quantum: kappa = h/m = {kappa:.1f} m^2/s   (BK at 0.6 eV: {kappa_bk:.0f})")
out = {}
for lam in (0.02, 0.035, 0.06):
    Om = lam*vf/(100*KPC)      # halo spin: v_rot ~ lambda v_f at ~100 kpc
    nv = 2*Om/kappa
    d = 1/np.sqrt(nv)
    Nv = nv*np.pi*RSF**2
    out[lam] = dict(Omega=Om, spacing_AU=d/1.496e11, N=Nv)
    print(f"  lambda={lam:.3f}: Omega={Om:.1e}/s  spacing={d/1.496e11:.3f} AU  N_vortices={Nv:.1e}")
lam = 0.035; Om = lam*vf/(100*KPC); nv = 2*Om/kappa; d = 1/np.sqrt(nv)
rho = 1e-21
eps = rho*kappa**2/(4*np.pi)*np.log(d/xi)          # line energy density
Etot = nv*np.pi*RSF**2 * eps * RSF                 # ~ column height R_SF
Ekin = 0.5*(1e12*MSUN)*cs**2
deficit = nv*np.pi*xi**2
print(f"energy: per-length {eps:.1e} J/m; total {Etot:.1e} J vs halo kinetic {Ekin:.1e} J -> ratio {Etot/Ekin:.1e}")
print(f"observability: fractional density deficit n_v*pi*xi^2 = {deficit:.1e}  (nothing)")
print(f"BK contrast at same Omega: spacing ratio ours/BK = sqrt(m_BK/m) = {np.sqrt(0.6/8.45):.2f} (x{np.sqrt(8.45/0.6):.1f} denser array)")
print("cores are NORMAL-PHASE FILAMENTS (X->0 on the axis): entropic stiffness off along threads of width ~xi")
print(f"  xi = {xi*1e3:.3f} mm; per-halo filament volume fraction ~ {deficit:.0e} — dynamically & observationally nil")

print("\n--- (2) the quintic vortex profile (extends [N11]) ---")
def vortex_profile(power):
    """shooting: u''+u'/r-u/r^2 = 2(u^power - u); u ~ a r at 0, u->1 at infinity."""
    def rhs(r, y): return [y[1], 2*(y[0]**power - y[0]) - y[1]/r + y[0]/r**2]
    def shoot(a):
        r0 = 1e-4
        sol = solve_ivp(rhs, [r0, 40], [a*r0, a], rtol=1e-10, atol=1e-12, dense_output=True, max_step=0.05)
        return sol
    lo, hi = 0.1, 2.0
    for _ in range(60):
        mid = 0.5*(lo+hi); s = shoot(mid); u_end = s.y[0][-1]
        if u_end > 1: hi = mid
        else: lo = mid
    s = shoot(0.5*(lo+hi))
    r = np.linspace(0.02, 12, 800); u = s.sol(r)[0]
    return r, np.clip(u, 0, 1.2), 0.5*(lo+hi)
rq, uq, aq = vortex_profile(5)
rc, uc, ac = vortex_profile(3)
# report radii in units of xi = hbar/(m c_s): equation units have xi_eq = 1/sqrt(2) (quintic), 1 (cubic)
r_half_q = rq[np.argmin(np.abs(uq - 1/np.sqrt(2)))] * np.sqrt(2)   # r/xi
r_half_c = rc[np.argmin(np.abs(uc - 1/np.sqrt(2)))]
print(f"quintic vortex: core slope a = {aq:.4f} (eq. units); density-half radius r_1/2 = {r_half_q:.3f} xi")
print(f"cubic   vortex: core slope a = {ac:.4f};              r_1/2 = {r_half_c:.3f} xi")
print(f"quintic core is {'narrower' if r_half_q < r_half_c else 'broader'} by {abs(1-r_half_q/r_half_c)*100:.0f}% at matched xi = hbar/(m c_s)")
print("  (same in-situ discriminant family as the soliton's broader-bottomed notch, Paper IV Fig. 1)")

json.dump(dict(kappa_m2_s=float(kappa), feynman=out,
               energy_ratio=float(Etot/Ekin), density_deficit=float(deficit),
               spacing_vs_BK="x%.1f denser at equal Omega" % np.sqrt(8.45/0.6),
               cores="normal-phase filaments, width ~ xi = %.3f mm" % (xi*1e3),
               quintic_vortex=dict(a=float(aq), r_half_xi=float(r_half_q)),
               cubic_vortex=dict(a=float(ac), r_half_xi=float(r_half_c))),
          open(os.path.join(DERIVED, "b2_vortices.json"), "w"), indent=1)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.plot(rq*np.sqrt(2), uq**2, "crimson", lw=1.7, label=f"quintic vortex ($r_{{1/2}}={r_half_q:.2f}\\,\\xi$)")
ax.plot(rc, uc**2, "#1f77b4", ls="--", lw=1.5, label=f"cubic vortex ($r_{{1/2}}={r_half_c:.2f}\\,\\xi$)")
ax.axhline(0.5, color=".6", lw=.6, ls=":")
ax.set(xlabel=r"$r/\xi$  ($\xi=\hbar/mc_s$)", ylabel=r"density $n/n_0$", xlim=(0, 8), ylim=(0, 1.05))
ax.legend(fontsize=7.5); ax.set_title("vortex core profiles at matched healing length", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures/fig_p4_vortex.png"), dpi=220, bbox_inches="tight")
print("\nwrote derived/b2_vortices.json, figures/fig_p4_vortex.png")