#!/usr/bin/env python
"""Paper I architecture figure: assumptions -> sectors -> the a0 lock -> phenomenology -> calendar."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(10.2, 6.2)); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
def box(x, y, w, h, text, fc, fs=8.3, ec=".25"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", fc=fc, ec=ec, lw=1))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs)
def arr(x1, y1, x2, y2, **kw):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=11,
                                 color=kw.get("c", ".35"), lw=kw.get("lw", 1.1)))
# inputs
box(0.4, 8.6, 4.2, 1.1, "[A2] Dorau–Much theorem (inherited)\n$S_{\\rm rel}$ = boost flux; $S=\\delta A/4$ ⇒ Einstein eqs,\n$\\alpha=8\\pi$, $\\Lambda$ an integration constant", "#dce9f7")
box(5.4, 8.6, 4.2, 1.1, "[A1] Berezhiani–Khoury architecture (inherited)\nnew boson, $m\\sim$ eV: superfluid in galaxies,\nnormal in clusters", "#dce9f7")
# postulates
box(0.4, 7.0, 4.2, 1.0, "[A3] fluctuation measure (this work)\n$P[\\varphi]\\propto e^{-S_{\\rm rel}[\\varphi]}$ — no free function", "#fdeacc")
box(5.4, 7.0, 4.2, 1.0, "[A4] acoustic entropy–area law (this work)\n$S_{\\rm rel}=\\delta A_{\\rm ac}/4G_{\\rm ph}$, $G_{\\rm ph}$ counted [N3]", "#fdeacc")
arr(2.5, 8.6, 2.5, 8.05); arr(7.5, 8.6, 7.5, 8.05)
# sectors
box(0.4, 4.9, 4.2, 1.6, "FUNDAMENTAL SECTOR  {A2, A3}\nvariance lemma $\\langle\\Delta K^2\\rangle{=}\\langle K\\rangle{\\Leftrightarrow}\\alpha{=}8\\pi$ [N1]\nseparation theorem, $\\sigma_\\Lambda/\\Lambda=1/\\sqrt{8\\pi}$ [N4]\nsecond-order law = Fisher info [N5]\n$\\Lambda$ = ledger constant: $w=-1$ identity", "#e8f2e2", fs=7.8)
box(5.4, 4.9, 4.2, 1.6, "EMERGENT SECTOR  {A1, A4}\nquintic class $\\mu\\propto n^2$: soliton, gap law,\nedge moment [N10–N14]\nquantum foundations: selection, conden-\nsation, resolution lock [N20–N22]", "#e8f2e2", fs=7.8)
arr(2.5, 7.0, 2.5, 6.55); arr(7.5, 7.0, 7.5, 6.55)
# the lock
box(2.6, 3.3, 4.8, 1.0, "ONE SHARED HORIZON (de Sitter)\n$a_0=c^2\\sqrt{\\Lambda/24\\pi}=(1.080\\pm0.010)\\times10^{-10}$ m s$^{-2}$  [N2]\nzero adjustable parameters; 4–11% below empirical", "#f6d9d9", fs=8.0)
arr(2.5, 4.9, 4.0, 4.35); arr(7.5, 4.9, 6.0, 4.35)
# phenomenology
box(0.4, 1.7, 9.2, 1.1, "PHENOMENOLOGY:  carrier ID $m=6$–$16$ eV, misalignment ALP [N7] · Bullet structure theorem [N8] · zero-parameter scorecard [N9]\n$\\sigma$-resolved RAR break measured, $m\\approx8$ eV [N16–N18] · derived tail (solar-safe) [N23] · adverse bar drag, recorded [N26]", "#eee6f5", fs=7.8)
arr(5.0, 3.3, 5.0, 2.85)
# calendar
box(0.4, 0.2, 9.2, 1.0, "FALSIFICATION CALENDAR [N15]:  Gaia DR4 wide binaries $\\delta\\tilde v=0.12$–$0.17$ — decided 2026-12-02 (kill if Newtonian)\nlensing = dynamics · $a_0(z)=$ const · $r<3\\times10^{-16}$ double kill · UV line 152–428 nm (positive channel)", "#fbf3d0", fs=7.8)
arr(5.0, 1.7, 5.0, 1.25)
fig.tight_layout()
fig.savefig("figures/fig_p1_architecture.pdf"); fig.savefig("figures/fig_p1_architecture.png", dpi=150)
print("wrote figures/fig_p1_architecture.pdf/.png")
