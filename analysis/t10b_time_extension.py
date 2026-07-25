#!/usr/bin/env python
"""
D29 addendum — time extension at the physical coupling (Theta ~ 0.85 target).

Question: is the marginal core coherence at 16 t_dyn an ENDPOINT or a snapshot
of slow entropic relaxation? (Unlike the bare particle, the entropic term is an
interaction: the sim HAS a condensation channel; its timescale is the unknown.)
Run 48000 steps (~48 t_dyn) at the gamma of the Theta-0.85 run; record the
core coherence C and granularity G every 4000 steps.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.ndimage import gaussian_filter
from t3_guard import DERIVED, ROOT

rng = np.random.default_rng(2910)                      # same seed: same collapse
NG = 256; L = 1.0; dx = L/NG
x = (np.arange(NG) - NG/2)*dx
X, Y = np.meshgrid(x, x, indexing="ij")
k1 = 2*np.pi*np.fft.fftfreq(NG, d=dx)
KX, KY = np.meshgrid(k1, k1, indexing="ij")
K2 = KX**2 + KY**2
K2inv = np.where(K2 > 0, 1/np.maximum(K2, 1e-30), 0.0)
A_GRAV = 6.7e5; R0 = 0.18; DT = 1.0e-6; NSTEP = 48000
KIN = np.exp(-0.5j*K2*DT)
GAMMA = json.load(open(os.path.join(DERIVED, "t10_collapse.json")))["theta_target_0.85"]["gamma"]

n0 = np.exp(-(X**2 + Y**2)/(2*R0**2)); n0 *= 1.0/(n0.sum()*dx*dx)
n0 *= (1 + 0.01*rng.standard_normal((NG, NG)))
psi = np.sqrt(np.abs(n0)).astype(np.complex128)

def poisson(n):
    return np.real(np.fft.ifft2(-np.fft.fft2(A_GRAV*(n - n.mean()))*K2inv))

def core_diag(psi):
    n = np.abs(psi)**2; m_tot = n.sum()
    cx = (X*n).sum()/m_tot; cy = (Y*n).sum()/m_tot
    r = np.hypot(X - cx, Y - cy)
    order = np.argsort(r.ravel()); csum = np.cumsum(n.ravel()[order])
    r_h = r.ravel()[order][np.searchsorted(csum, 0.5*csum[-1])]
    Ek = (0.5*np.abs(np.fft.ifft2(1j*KX*np.fft.fft2(psi)))**2
          + 0.5*np.abs(np.fft.ifft2(1j*KY*np.fft.fft2(psi)))**2).sum()*dx*dx
    v2 = 2*Ek/(m_tot*dx*dx); lam = 2*np.pi/np.sqrt(v2)
    n_s = gaussian_filter(n, float(np.clip(3*lam/dx, 3.0, r_h/(3*dx))))
    sel = r < r_h
    G = float((((n[sel] - n_s[sel])**2).mean())/max((n_s[sel]**2).mean(), 1e-30))
    C = float(np.abs(psi[sel].sum())**2/max((np.sqrt(n[sel])).sum()**2, 1e-30))
    nbar = float((n[sel]**2).sum()/n[sel].sum())
    return C, G, float((5/3)*2*GAMMA*nbar**2/v2)

print(f"=== D29 addendum: 48 t_dyn at gamma = {GAMMA:.3e} (the Theta~0.85 run) ===")
print(f"{'t/t_dyn':>8} {'C_core':>8} {'G_core':>8} {'Theta':>7}")
series = []
for s in range(NSTEP):
    n = np.abs(psi)**2
    psi *= np.exp(-0.5j*(poisson(n) + GAMMA*n**2)*DT)
    psi = np.fft.ifft2(KIN*np.fft.fft2(psi))
    n = np.abs(psi)**2
    psi *= np.exp(-0.5j*(poisson(n) + GAMMA*n**2)*DT)
    if (s+1) % 4000 == 0:
        C, G, th = core_diag(psi)
        series.append(dict(step=s+1, t_tdyn=round((s+1)*DT/1e-3, 1), C=round(C, 3), G=round(G, 3), theta=round(th, 2)))
        print(f"{(s+1)*DT/1e-3:8.0f} {C:8.3f} {G:8.3f} {th:7.2f}")

dC = series[-1]["C"] - series[3]["C"]
trend = "RISING (relaxation ongoing: marginal state is a snapshot, not an endpoint)" if dC > 0.1 else \
        ("FALLING" if dC < -0.1 else "FLAT (the marginal state is the endpoint at this Theta)")
print(f"\ncoherence trend after virialization: dC = {dC:+.3f}  ->  {trend}")
json.dump(dict(gamma=GAMMA, series=series, dC_late=float(dC), trend=trend),
          open(os.path.join(DERIVED, "t10b_time_extension.json"), "w"), indent=1)
print("wrote derived/t10b_time_extension.json")
