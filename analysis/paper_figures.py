#!/usr/bin/env python
"""
Publication figures for the five-manuscript portfolio (v0.2, hardened numbers).
Outputs PNG+PDF into figures/. Data sources: derived/*.json (frozen), plus the
in-script BdG snake computation (validated machinery, p4_referee_hardening.py).
  Letter : fig_letter_epsilon      (2-panel: eps vs H0 with bands; pure-Lambda)
  P2     : fig_p2_uncertainty_budget
  P3     : fig_p3_df2, fig_p3_phase_diagram, fig_p3_wb_dryrun, fig_p3_calendar
  P4     : fig_p4_profile, fig_p4_snake, fig_p4_edge_landscape, fig_p4_sequence
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2_guard import DERIVED, ROOT

FIG = os.path.join(ROOT, "figures"); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "legend.fontsize": 7.5, "figure.dpi": 110})
C = 2.99792458e8; MPC = 3.0856775814913673e22
def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=220, bbox_inches="tight")
    plt.close(fig); print(f"  wrote figures/{name}.png/.pdf")

# ================= Letter Fig 1 =================
def letter_epsilon():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 3.1))
    H = np.linspace(66.5, 74.6, 300); cH = C * H * 1e3 / MPC
    for a0, st, sy, c, lab in ((1.128e-10, 0.019e-10, 0.226e-10, "#1f77b4", "marginalized 1.128"),
                               (1.20e-10, 0.02e-10, 0.24e-10, "#444444", "canonical 1.20")):
        ax.plot(H, a0/cH, c=c, lw=1.6, ls="-" if c != "#444444" else "--", label=lab)
        ax.fill_between(H, (a0-st)/cH, (a0+st)/cH, color=c, alpha=.28, lw=0)
        ax.fill_between(H, (a0-sy)/cH, (a0+sy)/cH, color=c, alpha=.07, lw=0)
    eps_p = 0.16509
    ax.axhline(eps_p, color="crimson", lw=1.6)
    ax.fill_between(H, eps_p*0.9947, eps_p*1.0053, color="crimson", alpha=.25, lw=0)
    ax.text(66.8, eps_p+.0012, r"$\sqrt{\Omega_\Lambda/8\pi}$ (this work)", color="crimson", fontsize=8)
    ax.axhline(1/6, color="k", ls=":", lw=1); ax.text(73.6, 1/6+.0008, "1/6", fontsize=8)
    ax.axhline(1/(2*np.pi), color="k", ls="-.", lw=1); ax.text(73.4, 1/(2*np.pi)-.003, r"$1/2\pi$", fontsize=8)
    for h0, dh, lab in ((67.36, 0.54, "Planck"), (73.04, 1.04, "SH0ES")):
        ax.axvspan(h0-dh, h0+dh, color="gray", alpha=.18, lw=0); ax.text(h0, .1885, lab, ha="center", fontsize=7.5)
    ax.plot([70.3], [eps_p], "o", ms=5, mfc="none", mec="crimson", mew=1.4)
    ax.annotate(r"$H_0=70.3$", (70.3, eps_p), (71.6, .1568), fontsize=7.5, color="crimson",
                arrowprops=dict(arrowstyle="-", lw=.7, color="crimson"))
    ax.set(xlabel=r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]", ylabel=r"$\varepsilon=a_0/cH_0$",
           xlim=(66.5, 74.6), ylim=(.148, .192)); ax.legend(loc="lower left", framealpha=.9)
    ax.set_title(r"(a) across the $H_0$ tension")

    rows = [("prediction  $c^2\\sqrt{\\Lambda/24\\pi}$", 1.080, [.010], "crimson", "*", 11),
            ("modular  $c^2\\sqrt{\\Lambda/3}/2\\pi$", 0.862, [0], "gray", "s", 6),
            ("canonical RAR $g_\\dagger$", 1.20, [.02, .24], "#444444", "o", 5.5),
            ("marginalized fit", 1.128, [.019, .226], "#1f77b4", "o", 5.5)]
    for i, (lab, v, errs, c, mk, ms) in enumerate(rows):
        y = len(rows)-1-i
        if len(errs) > 1: bx.errorbar([v], [y], xerr=[[errs[1]], [errs[1]]], fmt="none", ecolor=c, alpha=.35, capsize=3)
        bx.errorbar([v], [y], xerr=[[errs[0]], [errs[0]]], fmt=mk, ms=ms, color=c, capsize=3, mew=1.2,
                    mfc="none" if mk == "s" else c)
        bx.text(1.47, y, lab, va="center", fontsize=8)
    bx.axvline(1.080, color="crimson", lw=.8, alpha=.5)
    bx.annotate("", (1.080, 2.55), (0.862, 2.55), arrowprops=dict(arrowstyle="<->", lw=.9))
    bx.text(0.968, 2.68, r"$\times\sqrt{\pi/2}$", ha="center", fontsize=8)
    bx.text(1.185, 0.28, "−4.4%", fontsize=7, color="#1f77b4"); bx.text(1.27, 1.28, "−11%", fontsize=7, color="#444")
    bx.set(xlabel=r"$a_0$  [$10^{-10}$ m s$^{-2}$]", xlim=(.78, 1.46), ylim=(-.6, 3.6), yticks=[])
    bx.set_title(r"(b) pure-$\Lambda$, zero parameters")
    fig.tight_layout(); save(fig, "fig_letter_epsilon")

# ================= P2 uncertainty budget =================
def p2_budget():
    d = json.load(open(os.path.join(DERIVED, "p2_hardening.json")))
    hs = [0.10, 0.15, 0.25]
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.axhspan(5.5, 6.2, color="crimson", alpha=.10, lw=0)
    ax.axhline(5.6, color="crimson", lw=1.3)
    ax.text(0.102, 5.66, r"quoted: $K_{95}=5.6\times10^{-13}$ (±12% envelope)", color="crimson", fontsize=7.5)
    for name, c, mk, lab in (("J1909-3744(PPTA)", "crimson", "o", "J1909−3744 (anchor, d err 1.2%)"),
                             ("J1640+2224(NG)", "gray", "s", "J1640+2224 (demoted: d err 31%)")):
        r = d[name]
        ax.plot(hs, [r[f"K95_uA|h={h}"]*1e13 for h in hs], mk+"-", color=c, ms=5, lw=1.2, label=lab+" — uniform-A")
        ax.plot(hs, [r[f"K95_lu|h={h}"]*1e13 for h in hs], mk+"--", color=c, ms=5, mfc="none", lw=1, alpha=.75,
                label=lab.split(" (")[0]+" — log-uniform")
    ax.set(xlabel=r"conditioning-kernel width $h$ (in $\gamma$)", ylabel=r"$K_{95}$  [$10^{-13}$ s$^2$ kpc$^{-1}$]",
           xticks=hs, ylim=(3.6, 9.9))
    ax.legend(loc="upper left", framealpha=.95)
    ax.set_title("ceiling uncertainty budget: anchors × kernel × prior")
    fig.tight_layout(); save(fig, "fig_p2_uncertainty_budget")

# ================= P3 figures =================
def p3_df2():
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.bar([0], [19.6], .55, color="crimson", alpha=.75, hatch="//", label="isolated MOND (excluded)")
    ax.bar([1], [9.3], .55, color="#1f77b4", alpha=.85, label="phonon-EFE (this work)")
    ax.errorbar([2], [8.5], yerr=[[3.1], [2.3]], fmt="ko", ms=6, capsize=4, label="observed (stellar $\\sigma$)")
    ax.set(xticks=[0, 1, 2], xticklabels=["isolated\nMOND", "phonon\nEFE", "observed"],
           ylabel=r"$\sigma$  [km s$^{-1}$]", ylim=(0, 23))
    ax.set_title("NGC1052-DF2 at $D=20$ Mpc")
    ax.text(.97, .03, "13 Mpc distance ⇒ test void (not failed)", transform=ax.transAxes,
            fontsize=7, color=".35", ha="right")
    ax.legend(loc="upper right", framealpha=.95)
    fig.tight_layout(); save(fig, "fig_p3_df2")

def p3_phase():
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    ax.set_xscale("log"); ax.set_xlim(4, 2500); ax.set_ylim(0, 5)
    ax.axvspan(4, 500, color="#1f77b4", alpha=.13, lw=0)
    ax.axvspan(500, 800, color="k", alpha=.10, lw=0, hatch="///")
    ax.axvspan(800, 2500, color="orange", alpha=.10, lw=0)
    ax.text(38, 4.45, "SUPERFLUID PHASE\nphonon force active — MOND phenomenology", ha="center", fontsize=8, color="#1f77b4")
    ax.text(632, 0.55, "$\\sigma_{\\rm crit}$\n$m=6$–$16$ eV", ha="center", fontsize=7.5)
    ax.text(1380, 4.45, "NORMAL PHASE\ncondensate = CDM\n(MOND dies)", ha="center", fontsize=7.5, color="darkorange")
    sysrows = [("dwarf spheroidals", 5, 30, 3.6), ("disc galaxies", 80, 300, 2.8),
               ("galaxy groups (sharpest test)", 250, 750, 2.0), ("clusters", 700, 1500, 1.2)]
    for lab, lo, hi, y in sysrows:
        ax.plot([lo, hi], [y, y], lw=5, color=".25", solid_capstyle="round", alpha=.75)
        ax.text(np.sqrt(lo*hi), y+.16, lab, ha="center", fontsize=7.5)
    ax.annotate(r"$T/T_c \propto M^{2/3}$: smaller systems are more superfluid",
                (5.5, .45), fontsize=7.5, color=".3")
    ax.set(xlabel=r"velocity dispersion $\sigma$  [km s$^{-1}$]", yticks=[])
    ax.set_title("one phase boundary: MOND-in-galaxies, CDM-in-clusters")
    fig.tight_layout(); save(fig, "fig_p3_phase_diagram")

def p3_wb():
    d0 = json.load(open(os.path.join(DERIVED, "t2_5_results.json")))["bins"]
    d = json.load(open(os.path.join(DERIVED, "t2_5b_repair.json")))["bins"]
    order = [("control_0.2-1kAU", 0.45), ("2-5kAU", 3.2), ("5-10kAU", 7.1), ("10-30kAU", 17.3)]
    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    ax.set_xscale("log")
    ax.axhspan(0.124, 0.167, xmin=.28, color="green", alpha=.18, lw=0)
    ax.text(28, .145, "frozen\nprediction\n$0.124$–$0.167$", fontsize=7.5, color="darkgreen", ha="center", va="center")
    ax.axvspan(.2, 1.0, color="#1f77b4", alpha=.06, lw=0)
    ax.axhline(0, color="k", lw=.7)
    for est, mk, c, off in (("thermal", "o", "#1f77b4", .93), ("hwang", "^", "#7f4fc9", 1.075)):
        x = [p*off for _, p in order]
        ax.plot(x, [d0[b][est]["delta_v"] for b, _ in order], mk, mfc="none", mec=".65", ms=4.5)
        y = [d[b][est]["delta_v"] for b, _ in order]
        e = [d[b][est]["delta_v_err"] for b, _ in order]
        ax.errorbar(x, y, yerr=e, fmt=mk, color=c, ms=5.5, capsize=3, lw=1, label=f"{est} e-prior (repaired)")
    ax.plot([], [], "o", mfc="none", mec=".65", ms=4.5, label="pre-repair dry-run")
    ax.annotate("control repaired (T2.5b): $-0.006\\pm0.004$\n$\\eta$, $f_t$ frozen for DR4",
                (.5, -.009), (1.15, .028), fontsize=7, color="#1f77b4", va="top",
                arrowprops=dict(arrowstyle="-", lw=.8, color="#1f77b4"))
    ax.set(xlabel=r"projected separation $r_p$  [kAU]", ylabel=r"$\delta\tilde v$",
           ylim=(-.04, .27), xlim=(.28, 40))
    ax.set_xticks([0.5, 1, 2, 5, 10, 30]); ax.set_xticklabels(["0.5", "1", "2", "5", "10", "30"])
    ax.legend(loc="upper left", framealpha=.95, fontsize=6.5)
    ax.set_title("eDR3 after control calibration — no adjudication; DR4 runs this code", fontsize=8.5)
    fig.tight_layout(); save(fig, "fig_p3_wb_dryrun")

def p3_calendar():
    fig, ax = plt.subplots(figsize=(6.8, 3.3))
    rows = [  # (label, start, end|None(point)|arrow, color, note)
        ("public deposit (protocol step)", 2026.6, 2026.9, ".4", "pre-registration force"),
        ("Gaia DR4 wide binaries", 2026.92, None, "crimson", "KILL if Newtonian — 2026-12-02"),
        ("$\\sigma$-resolved RAR / groups", 2027.0, 2030.5, "darkorange", "break must track $\\sigma$"),
        ("lensing = dynamics", 2026.6, 2032.4, "darkorange", "standing falsifier"),
        ("cluster collisions $\\sigma/m$", 2026.6, 2032.4, "darkorange", ""),
        ("CMB-S4 tensor $r$", 2029.0, 2032.4, "crimson", "any detection = DOUBLE KILL"),
        ("UV decay line 152–428 nm", 2027.0, 2032.4, "green", "positive channel"),
        ("cold-atom $C_T$ (Paper IV)", 2027.5, 2032.4, "green", "positive channel"),
    ]
    for i, (lab, s, e, c, note) in enumerate(rows):
        y = len(rows)-i
        if e is None:
            ax.plot([s], [y], "D", color=c, ms=8)
        else:
            ax.plot([s, e], [y, y], lw=5, color=c, alpha=.55, solid_capstyle="round")
            if e >= 2032.3: ax.annotate("", (e+.25, y), (e-.02, y), arrowprops=dict(arrowstyle="->", color=c, alpha=.7))
        ax.text(2026.5, y, lab+"  ", ha="right", va="center", fontsize=8)
        if note: ax.text((s if e is None else s)+.12, y+.33, note, fontsize=6.7, color=c)
    ax.axvline(2026.92, color="crimson", lw=.8, ls=":", alpha=.6)
    ax.set(xlim=(2024.9, 2032.9), ylim=(.3, len(rows)+1.0), yticks=[],
           xticks=[2026, 2027, 2028, 2029, 2030, 2031, 2032])
    ax.set_title("the falsification calendar — kill tests (red), consistency (amber), discovery (green)")
    for sp in ("left", "right", "top"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); save(fig, "fig_p3_calendar")

# ================= P4 figures =================
def quintic_u(z):  # z in units of xi
    t = np.tanh(z); return np.sqrt(2)*t/np.sqrt(3-t**2)

def p4_profile():
    z = np.linspace(-5, 5, 800)
    fig, ax = plt.subplots(figsize=(4.6, 3.1))
    ax.plot(z, quintic_u(z), "crimson", lw=1.8, label=r"quintic: $\sqrt{2}\tanh/\sqrt{3-\tanh^2}$")
    ax.plot(z, np.tanh(z), "#1f77b4", ls="--", lw=1.5, label=r"cubic: $\tanh(z/\xi)$")
    ax.axhline(0, color="k", lw=.5)
    ax.text(-4.75, .28, "nodal slopes:\nquintic $\\sqrt{2/3}$, cubic $1$\n(quintic broader-bottomed)", fontsize=7.3, color=".25")
    ins = ax.inset_axes([.60, .12, .36, .34])
    ins.plot(z, quintic_u(z)**2, "crimson", lw=1.4); ins.plot(z, np.tanh(z)**2, "#1f77b4", ls="--", lw=1.1)
    ins.set(xlim=(-3, 3), ylim=(0, 1.05)); ins.set_title(r"density $n/n_0$", fontsize=7); ins.tick_params(labelsize=6)
    ax.set(xlabel=r"$z/\xi$", ylabel=r"$\psi/\sqrt{n_0}$", xlim=(-5, 5), ylim=(-1.15, 1.15))
    ax.legend(loc="upper left", framealpha=.95)
    ax.set_title("the black soliton at matched healing length")
    fig.tight_layout(); save(fig, "fig_p4_profile")

def p4_snake():
    def scan(kind, N=600, L=44.0):
        zg = np.linspace(-L/2, L/2, N); h = zg[1]-zg[0]
        D2 = (np.diag(np.full(N-1, 1.), -1) - 2*np.eye(N) + np.diag(np.full(N-1, 1.), 1))/h**2
        if kind == "cubic":
            psi = np.tanh(zg); xiph = 1.0; Vp, Vm = 3*psi**2, psi**2
        else:
            tt = np.tanh(zg/(1/np.sqrt(2))); psi = np.sqrt(2)*tt/np.sqrt(3-tt**2)
            xiph = 1/np.sqrt(2); Vp, Vm = 5*psi**4, psi**4
        ks = np.arange(0.03, 1.16, 0.02)/xiph
        gam = []
        for k in ks:
            Lp = -0.5*(D2 - k**2*np.eye(N)) - np.eye(N) + np.diag(Vp)
            Lm = -0.5*(D2 - k**2*np.eye(N)) - np.eye(N) + np.diag(Vm)
            w2 = np.linalg.eigvals(Lm @ Lp)
            gam.append(np.sqrt(max(0., -np.min(w2.real))))
        return ks*xiph, np.array(gam)
    fig, ax = plt.subplots(figsize=(4.8, 3.2))
    for kind, c, w in (("quintic", "crimson", .751), ("cubic", "#1f77b4", 1.000)):
        kx, g = scan(kind)
        ax.plot(kx, g, color=c, lw=1.8, ls="-" if kind == "quintic" else "--", label=kind)
        ax.axvline(w, color=c, lw=.8, ls=":", alpha=.8)
        ax.text(w, .335, f"{w:.3f}", color=c, ha="center", fontsize=7.5)
    ax.axhline(0.3, color=".4", lw=.8, ls="-."); ax.text(.045, .306, "pinning threshold 0.3 $\\mu/\\hbar$", fontsize=7, color=".35")
    ax.annotate(r"$0.302$ @ $k\xi=0.48$", (.48, .302), (.12, .27), fontsize=7.5, color="crimson",
                arrowprops=dict(arrowstyle="-", lw=.7, color="crimson"))
    ax.set(xlabel=r"$k_\parallel \xi$  ($\xi=\hbar/mc_s$)", ylabel=r"snake growth  ${\rm Im}\,\omega$  [$\mu/\hbar$]",
           xlim=(0, 1.15), ylim=(0, .36))
    ax.legend(loc="upper right", framealpha=.95)
    ax.set_title("transverse instability: quintic vs cubic (code validated at 1.000)")
    fig.tight_layout(); save(fig, "fig_p4_snake")

def p4_edge():
    m = json.load(open(os.path.join(DERIVED, "b2_mode_selection.json")))
    tgt = m["target"]["k1_over_n13"]
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.axvline(tgt, color="crimson", lw=1.5)
    ax.axvspan(tgt*.98, tgt*1.02, color="crimson", alpha=.12, lw=0)
    ax.text(tgt, 8.35, "required 4.257\n(±2% measurement goal)", ha="center", fontsize=7.5, color="crimson")
    entries = [("sharp edge + rms", m["movement3"]["edges"]["sharp edge, rms selection"]["k1"], ".5"),
               ("Debye $\\nu=1$ (= sharp Debye edge)", m["movement2"]["candidates"]["nu=1 (Debye)"]["k1"], ".35"),
               ("Gaussian edge + mean", m["movement3"]["edges"]["Gaussian edge, mean selection"]["k1"], ".5"),
               ("$S(k)$-weighted budget (drifts)", None, "#7f4fc9"),
               ("$\\nu=4/3$ (numerological)", m["movement2"]["candidates"]["nu=4/3"]["k1"], ".5"),
               ("Gaussian edge + rms  (best, +2.0%)", m["movement3"]["edges"]["Gaussian edge, rms selection"]["k1"], "darkgreen"),
               ("$\\nu=3/2$", m["movement2"]["candidates"]["nu=3/2"]["k1"], ".6"),
               ("$\\nu=2$", m["movement2"]["candidates"]["nu=2 (two chiralities)"]["k1"], ".6")]
    for i, (lab, v, c) in enumerate(entries):
        y = len(entries)-i
        if v is None:
            g, co = m["movement1"]["k1_over_n13_by_env"]["galactic"], m["movement1"]["k1_over_n13_by_env"]["cosmic"]
            ax.annotate("", (co, y), (g, y), arrowprops=dict(arrowstyle="->", color=c, lw=1.6))
            ax.text((g+co)/2, y+.30, "3.97→3.90", fontsize=6.6, color=c, ha="center")
        else:
            ax.plot([v], [y], "o", ms=6.5, color=c)
            ax.text(v, y+.30, f"{(v/tgt-1)*100:+.1f}%", ha="center", fontsize=6.6, color=c)
        ax.text(2.86, y, lab, ha="right", va="center", fontsize=7.6)
    ax.set(xlim=(2.9, 5.15), ylim=(.3, 9.4), yticks=[],
           xlabel=r"branch-edge budget  $k_1/n^{1/3}$   (deviations in $k$; in $C_T$ ≈ 2×)")
    ax.set_title("no counting closes the edge — the moment is a material constant", fontsize=9)
    for sp in ("left", "right", "top"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); save(fig, "fig_p4_edge_landscape")

def p4_sequence():
    fig, axs = plt.subplots(1, 3, figsize=(8.0, 2.7))
    x = np.linspace(-6, 6, 500)
    for k, (ax, ttl) in enumerate(zip(axs, ["(i) imprint planar soliton\n(box trap)",
                                            "(ii) pin: blue-detuned sheet\nbarrier ≳ 0.3 μ at the node",
                                            "(iii) Bragg: branch → edge\nreport ⟨k²⟩$^{1/2}_E/k_g$"])):
        ax.set_title(ttl, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
        if k < 2:
            ax.plot(x, quintic_u(x)**2, "crimson", lw=1.6)
            ax.add_patch(plt.Rectangle((-6, 0), 12, 1.25, fill=False, ec=".3", lw=1.2))
            ax.set(xlim=(-6.6, 6.6), ylim=(-.12, 1.55))
            ax.text(0, -.085, "n(z)", ha="center", fontsize=7, color="crimson")
            if k == 1:
                ax.plot(x, 1.28*np.exp(-x**2/.09), color="#1f77b4", lw=1.6)
                ax.text(0, 1.40, "light sheet", ha="center", fontsize=7, color="#1f77b4")
                ax.annotate("", (0.85, .55), (2.3, .55), arrowprops=dict(arrowstyle="->", color=".4"))
                ax.annotate("", (-0.85, .55), (-2.3, .55), arrowprops=dict(arrowstyle="->", color=".4"))
                ax.text(3.15, .52, "snake\nsuppressed", fontsize=6.4, color=".4", va="center")
        else:
            kk = np.linspace(0, 1.6, 300)
            om = kk*np.sqrt(1+kk**2/2)
            kg = 1.05
            for a, b in zip(kk[:-1], kk[1:]):
                al = 1.0 if a < kg*.8 else max(0., 1-((a-kg*.8)/(.55*kg)))
                ax.plot([a, b], [a*np.sqrt(1+a**2/2), b*np.sqrt(1+b**2/2)], color="crimson", lw=1.7, alpha=al)
            ax.fill_between(kk, om*1.15+.25, 2.9, where=kk > .55, color=".8", alpha=.5, lw=0)
            ax.text(1.18, 2.55, "multi-excitation\ncontinuum", fontsize=6.4, color=".4", ha="center")
            ax.axvline(kg, color=".4", lw=.7, ls=":"); ax.text(kg, -.28, r"$k_g$", ha="center", fontsize=7)
            ke = np.linspace(.7, 1.5, 200)
            ax.plot(ke, 0.75*np.exp(-(ke-1.26*kg)**2/(2*.14**2)), color="darkgreen", lw=1.4)
            ax.text(1.32, .84, r"$E(k)=-dZ/dk$", fontsize=6.6, color="darkgreen")
            ax.annotate("", (1.26*kg, .1), (1.26*kg, .55), arrowprops=dict(arrowstyle="->", color="darkgreen", lw=1))
            ax.text(1.26*kg+.03, .16, r"$1.20\,k_g$?", fontsize=6.6, color="darkgreen")
            ax.set(xlim=(0, 1.65), ylim=(-.35, 2.9)); ax.set_xlabel("k", fontsize=7); ax.set_ylabel(r"$\omega$", fontsize=7)
    fig.tight_layout(); save(fig, "fig_p4_sequence")

print("=== generating portfolio figures ===")
letter_epsilon(); p2_budget(); p3_df2(); p3_phase(); p3_wb(); p3_calendar()
p4_profile(); p4_snake(); p4_edge(); p4_sequence()
print("done.")