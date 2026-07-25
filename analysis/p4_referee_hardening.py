#!/usr/bin/env python
"""
Paper IV referee-hardening:
 (1) closed-form soliton: direct ODE residual, central slope sqrt(2/3), decay rate 2/xi
 (2) BdG CODE VALIDATION on the cubic black soliton (known snake window ~1 in
     xi = hbar/m c_s units) run with the SAME machinery that produced the quintic
     0.751 window -- the referee's "did you validate on the known case?"
 (3) T_c-as-density-scale identity: k_B T_c/hbar = (2 pi hbar/m)(n/zeta(3/2))^(2/3)
     -- C_T notation invokes NO thermal physics (kills the thermal-gap misreading)
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2_guard import DERIVED
from scipy.special import zeta

print("=== (1) closed-form soliton checks (units: hbar=m=1, gamma=1, n0=1, mu=1) ===")
cs = np.sqrt(2.0); xi = 1 / cs
z = np.linspace(-14, 14, 20001); dz = z[1] - z[0]
t = np.tanh(z / xi)
u = np.sqrt(2) * t / np.sqrt(3 - t**2)
upp = np.gradient(np.gradient(u, dz), dz)
res = -0.5 * upp + u**5 - u
i = slice(2000, 18001)
print(f"max |ODE residual| (interior) = {np.max(np.abs(res[i])):.1e}")
slope = np.gradient(u, dz)[len(z)//2]
print(f"central slope = {slope:.6f}  vs sqrt(2/3)/xi = {np.sqrt(2/3)/xi:.6f}")
j = (z > 3) & (z < 9)
decay = -np.polyfit(z[j], np.log(1 - u[j]), 1)[0]
print(f"asymptotic decay rate = {decay:.4f}  vs 2/xi = {2/xi:.4f}")

print("\n=== (2) BdG snake windows: cubic validation + quintic, same code ===")
def snake(kind, N=700, L=46.0):
    zg = np.linspace(-L/2, L/2, N); h = zg[1] - zg[0]
    D2 = (np.diag(np.full(N-1, 1.0), -1) - 2*np.eye(N) + np.diag(np.full(N-1, 1.0), 1)) / h**2
    if kind == "cubic":     # psi = tanh(z), mu = 1, c_s = 1, xi_phys = 1
        psi = np.tanh(zg); mu = 1.0; xiph = 1.0
        Vp, Vm = 3*psi**2, psi**2
    else:                   # quintic: mu = 1, c_s = sqrt(2), xi_phys = 1/sqrt(2)
        tt = np.tanh(zg / (1/np.sqrt(2)))
        psi = np.sqrt(2)*tt/np.sqrt(3 - tt**2); mu = 1.0; xiph = 1/np.sqrt(2)
        Vp, Vm = 5*psi**4, psi**4
    out = []
    for k in np.arange(0.02, 1.9, 0.04):
        Lp = -0.5*(D2 - k**2*np.eye(N)) - mu*np.eye(N) + np.diag(Vp)
        Lm = -0.5*(D2 - k**2*np.eye(N)) - mu*np.eye(N) + np.diag(Vm)
        w2 = np.linalg.eigvals(Lm @ Lp)
        gmax = np.sqrt(max(0.0, -np.min(w2.real)))
        out.append((k*xiph, gmax))
    out = np.array(out)
    un = out[out[:, 1] > 1e-3]
    kedge = un[-1, 0] + (out[out[:, 0] > un[-1, 0]][0, 0] - un[-1, 0])/2 if len(un) else 0
    ipk = np.argmax(out[:, 1])
    return kedge, out[ipk, 1], out[ipk, 0]
for kind, known in (("cubic", "literature ~1 (validation)"), ("quintic", "D20: 0.751, 0.302 @ 0.48")):
    ke, gm, kp = snake(kind)
    print(f"{kind:8s}: window k*xi < {ke:.3f}, max growth {gm:.3f} mu/hbar at k*xi = {kp:.3f}   [{known}]")

print("\n=== (3) T_c as a pure density scale ===")
print("k_B T_c/hbar = (2 pi hbar/m) (n/zeta(3/2))^(2/3); check dimensionless ratio at two densities:")
for n in (1e12, 1e20):     # SI-ish arbitrary
    hbar_m = 1.0           # units out; ratio test
    tc = 2*np.pi * (n/zeta(1.5, 1))**(2/3)   # zeta(3/2)=2.612
    print(f"  n = {n:.0e}: (k_B T_c/hbar)/n^(2/3) = {tc/n**(2/3):.6f}  (constant = 2 pi/zeta(3/2)^(2/3))")
print(f"  2 pi / zeta(3/2)^(2/3) = {2*np.pi/2.6123753486854883**(2/3):.6f} -- no thermal physics, pure notation")

json.dump(dict(ode_residual=float(np.max(np.abs(res[i]))), slope=float(slope),
               decay=float(decay)), open(os.path.join(DERIVED, "p4_hardening.json"), "w"), indent=1)
print("\nsaved derived/p4_hardening.json")