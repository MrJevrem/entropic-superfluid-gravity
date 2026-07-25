# Entropic Gravity on a Dark Superfluid — analysis code and frozen prediction ledgers

Code, derived products, and pre-registered prediction ledgers for a six-manuscript series deriving the Milgrom acceleration scale a₀ = c²√(Λ/24π) from the Dorau–Much relative-entropy theorem (PRL **136**, 091602) joined to Berezhiani–Khoury superfluid dark matter (PRD **92**, 103510) by two postulates — and staking the result on a **frozen, single-fork Gaia DR4 wide-binary prediction (δṽ = 0.12–0.17) decided on 2026-12-02**.

Manuscripts: *(links added at deposit — Zenodo DOI / arXiv IDs)*

| Paper | Content |
|---|---|
| I (flagship) | assumptions [A1–A4], the a₀ law, counted G_ph, dark-energy sector, master results ledger [N1–N19] |
| II | the pulsar-timing ceiling K₉₅ = 5.6×10⁻¹³ s² kpc⁻¹ (self-contained archival measurement) |
| III | dark-matter phenomenology, carrier identity card, falsification calendar |
| IV | the quintic superfluid: exact black soliton, branch-edge physics, laboratory analog |
| V | the σ-resolved radial-acceleration relation: the break, its baryonic ceiling, the carrier-mass band |
| Letter | the a₀ = c²√(Λ/24π) result in PRL format |

## What is (and is not) in this repository

- `analysis/` — every analysis script (guards, freezes, pipelines, derivation chain, hardening suites, figures).
- `derived/` — frozen prediction ledgers (`*_locked_predictions*.json`) and per-step JSON results. **The ledgers are the pre-registration record**: each embeds the UTC timestamp and development-repository commit hash at freeze time, and revisions only ever append versioned entries with logged reasons.
- **Not included:** archival datasets (large; all public — download table below), manuscript sources, and internal working documents. Figures regenerate from `analysis/paper_figures.py`.

## Data: what is used, where to get it

Scripts expect the layout below (paths relative to the repository root). Nothing here is proprietary; every dataset is a published archive.

| Dataset | Used by | Source | Local path |
|---|---|---|---|
| NANOGrav 15 yr narrowband v2.1.0 (single-pulsar noise chains, par files) | T1.1/Paper II | NANOGrav data release (Zenodo) | `data/T1.1_pulsar_timing/nanograv_15yr/` |
| NANOGrav 15 yr stochastic-analysis cores (joint CURN chain) | T1.1 v4 | `nanograv/15yr_stochastic_analysis` (GitHub, presampled cores; commit-pinned) | same tree |
| NANOGrav 15 yr KDE free spectra v1.1.0 | T1.1 robustness | NANOGrav release (Zenodo) | same tree |
| PPTA DR3 joint CURN chain (`chain_commonNoise_pl_nocorr_freegam_DE440.npy`) | T1.1 v3/v4 | PPTA DR3 analysis repository (commit-pinned) | `data/T1.1_pulsar_timing/ppta_dr3/` |
| EPTA DR2 maximum-likelihood noise solutions | T1.1 v3 | EPTA DR2 release | `data/T1.1_pulsar_timing/epta_dr2/` |
| ATNF pulsar catalogue v2.8.1 | distances | ATNF psrcat | `data/T1.1_pulsar_timing/` |
| SPARC mass models | T2/T3 | SPARC database (Lelli, McGaugh & Schombert) | `data/T2/sparc/` |
| El-Badry+21 eDR3 wide-binary catalog (`all_columns_catalog.fits.gz`, ~1.3 GB) | T2.5/T2.5b | catalog release (Zenodo) | `data/T2/wb/` |
| Sun+09 Chandra group tables | T3 | ApJ **693**, 1142 e-print/CDS tables | `data/T3/sun09/` |
| Vikhlinin+06 relaxed-cluster fits | T3 | astro-ph/0507092 e-print source | `data/T3/vikhlinin06/` |
| ACCEPT cluster tables (`table1.dat.gz`, `table5.dat.gz`) | T3 | ACCEPT archive | `data/T3/accept/` |
| SLUGGS/Alabi+16 tracer-mass tables | T3 | MNRAS **460**, 3838 e-print tables | `data/T3/sluggs/` |

Planck 2018 values enter as locked constants inside the scripts (no download).

**Staging rule (protocol-relevant):** the guards treat file modification times as staging times. Stage data by direct download *after* the relevant freeze; if extracting archives that restore historical mtimes, `touch` the extracted tree and record the download event (see the guard docstrings).

## The pre-registration protocol

`t2_guard.py` / `t3_guard.py` are imported by every comparison script and refuse to run unless: the frozen ledger exists, the analysis tree is git-clean, the freeze commit is an ancestor of HEAD, and every comparison-data file was staged *after* the freeze. Freeze scripts (`t2_predictions_freeze.py`, `t3_predictions_freeze.py`) compute and lock predictions — including both branches wherever a derivation forks — before any comparison data is staged. Revisions are versioned, reasons logged, superseded versions never deleted (see the `revisions` blocks inside the ledgers, including the pre-data wide-binary coefficient correction and the T2.5b control-repair nuisances frozen for Gaia DR4).

This public snapshot carries the ledgers and their embedded timestamps; the commit hashes inside them refer to the private development history, whose full ancestry is available for inspection on request. The externally binding timestamps are the public deposit dates (DOI / arXiv).

**Frozen headline predictions:** Gaia DR4 wide binaries δṽ = 0.12–0.17 (single fork; kill if Newtonian; 2026-12-02) · PTA ceiling K₉₅ = 5.6×10⁻¹³ s² kpc⁻¹ · condensation break σ_crit(m) window [250, 1000] km s⁻¹ with the measured pin m ≈ 8 eV · primordial tensors r < 3×10⁻¹⁶ (any detection is a double kill).

## Script guide

**Guards and freezes** — `t2_guard.py`, `t3_guard.py` (enforcement); `t2_predictions_freeze.py`, `t3_predictions_freeze.py` (the freezes).

**T1.1 — pulsar-timing search and ceiling (Paper II):** `t1_1_pulsar_distance_regression.py` (the naive estimator and its apparent signal), `t1_1_v2_common_gamma_envelope.py` (audit-corrected: index conditioning + envelope), `t1_1_v3_cross_pta.py` (NANOGrav × PPTA DR3 × EPTA DR2), `t1_1_v4_joint_vs_joint.py` (method-matched consistency; the 12/25 reabsorption), `t1_1_freespectrum_robustness.py`, `t1_1_method_comparison_figure.py`.

**T2 — branch discrimination (Papers I/III):** `t2_0_squeeze.py` (decoherence/squeeze trade-off atlas), `t2_1_compare_ceiling.py` (frozen branch predictions vs the ceiling), `t2_2_sightline_entropy.py`, `t2_3_early_universe.py` (ΔN_eff, FIRAS gates), `t2_4_a0_of_z.py` (the a₀(z) lock), `t2_5_dr3_dryrun.py` (eDR3 wide-binary dry-run, both contested estimator families), `t2_5b_control_repair.py` (the control-bin repair; nuisances frozen for DR4), `t2_survival.py` (assembled verdicts).

**T3 — the σ-resolved RAR (Paper V):** `t3_1_anchors.py` through `t3_5_systematics.py` (anchors, shelf, break fit, discriminator grid, fork forest), `t3_6_fb_baseline.py` (the partial-restoration f_b(σ) model and the m = 7.9 eV pin), `t3_8_synthesis.py` (verdicts vs the frozen ledger).

**D-chain — conditional derivations (the theory record, D0–D26):** `b2_derivations.py` (the base chain), `b2_epsilon_decision.py` (1/6 vs 1/2π), `b2_gph_counting.py` (G_ph), `b2_lambda_bookkeeping.py` / `b2_f_floor_test.py` / `b2_ossw_saturation.py` (the separation-theorem exhaustion), `b2_theta2_identity.py` / `b2_consistency_triad.py` (the second-order law), `b2_dm_implications.py` / `b2_production.py` (carrier identity), `b2_x32_coexistence.py`, `b2_cw_matching.py`, `b2_screen_eigenvalue.py` / `b2_notch_bdg.py` / `b2_mode_selection.py` (the C_T program), `b2_phonon_stability.py` (D22), `b2_vortices.py` (D23), `b2_rsf_outer_halo.py` (D24), `b2_isolated_giants.py` (D25), `b2_quantum_breaking.py` (D26).

**Hardening and figures:** `p1_referee_hardening.py` (Paper I checks + the portfolio lint), `p2_…`–`p4_…`, `p_letter_hardening.py` (per-paper claim verification), `paper_figures.py` (all publication figures).

## Reproducing

Python ≥ 3.11 with `numpy`, `scipy`, `astropy`, `matplotlib`, `pandas` (T1.1 chain handling additionally uses `la_forge` and `ceffyl`-format products). Each script is standalone (`python analysis/<script>.py`) and writes its JSON record to `derived/`; comparison scripts self-verify the protocol via the guards. `analysis/paper_figures.py` regenerates every figure.

## AI contribution

This research program was carried out with the help of **Claude Fable 5**, an AI research system built by Anthropic, which performed the derivations, wrote and executed the analysis code, and drafted the manuscripts; the human author (Marko Jevremovic) directed the program and takes sole responsibility for the content. Current publication policies do not permit AI systems to hold authorship; each manuscript carries a contribution statement recording that co-authorship was intended.

## License and citation

Code: MIT (see `LICENSE`). Citation metadata: `CITATION.cff` (the Zenodo concept DOI is added there at first release). Public repository: <https://github.com/MrJevrem/entropic-superfluid-gravity> — a snapshot of `analysis/` + `derived/` mirrored by `sync_public.sh` from the private development repository.
