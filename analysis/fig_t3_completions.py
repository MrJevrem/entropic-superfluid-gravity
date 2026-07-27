#!/usr/bin/env python
"""Paper V fig: (a) the budget-completion ladder for sigma_b; (b) the cluster amplitude by method."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.4))
# (a) sigma_b ladder
rows = [("baseline (BCG stars,\nmeasured gas)", 254, 215, 315, ".35"),
        ("stellar completion\n(T3.9, bg-subtracted)", 508, None, None, "#7f4fc9"),
        ("gas restoration\n(T3.6, retentive baseline)", 657, 621, 704, "#1f77b4"),
        ("combined (T3.10)", 628, 599, 679, "crimson"),
        ("combined + LF\ncompletion", 638, 608, 697, "crimson")]
a.axvspan(476, 756, color="orange", alpha=.15, lw=0)
a.axvline(600, color="darkorange", ls=":", lw=1.4)
a.text(600, 4.55, " frozen fiducial 600\n (m = 8.45 eV; ρ band 476–756)", fontsize=7, color="darkorange")
for i, (lab, v, lo, hi, c) in enumerate(rows):
    xerr = None if lo is None else [[v-lo], [hi-v]]
    a.errorbar([v], [i], xerr=xerr, fmt="o", color=c, ms=6, capsize=4, lw=1.4)
    a.text(150, i, lab, fontsize=7.2, va="center", ha="left")
a.set(xlim=(140, 900), ylim=(-.6, 5.1), yticks=[], xlabel=r"break location $\sigma_b$ [km s$^{-1}$]")
a.set_title("(a) budget completions walk the break onto the fiducial", fontsize=8.6)
# (b) amplitude by method
mrows = [("X-ray chain (V06)", .44, None, "#d62728"),
         ("X-ray, cosmic $f_b$\nconcession", .374, None, "#d62728"),
         ("caustic ⊗ X-ray,\nσ>1000 (T3.7)", .216, .041, "#1f77b4"),
         ("caustic, SZ-mass fork", .175, None, "#1f77b4"),
         ("caustic, c=5 fork", .147, None, "#1f77b4"),
         ("caustic ⊗ X-ray, all\n(Duffy c)", .113, .032, "#1f77b4")]
b.axvline(.15, color=".4", ls="--", lw=1)
b.text(.152, 5.0, "frozen P2\nthreshold", fontsize=7, color=".35")
for i, (lab, v, e, c) in enumerate(mrows):
    b.errorbar([v], [i], xerr=None if e is None else [[e], [e]], fmt="s", color=c, ms=5.5, capsize=4, lw=1.2)
    b.text(-.02, i, lab, fontsize=7.0, va="center", ha="right")
b.set(xlim=(-.24, .52), ylim=(-.6, 5.6), yticks=[], xlabel=r"cluster deviation $D$ [dex]")
b.set_title("(b) the normal-phase amplitude is method-dependent", fontsize=8.6)
fig.tight_layout()
fig.savefig("figures/fig_t3_completions.png", dpi=200, bbox_inches="tight")
fig.savefig("figures/fig_t3_completions.pdf", bbox_inches="tight")
print("wrote figures/fig_t3_completions.png/.pdf")
