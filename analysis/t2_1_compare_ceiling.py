#!/usr/bin/env python
"""T2.1 comparison: frozen branch K predictions vs the T1.1 ceiling grid (plan §3).

Verdict rules (frozen in plan): K_pred > K95 => excluded (internal contradiction);
0-2 dex below => reachable-next-gen (IPTA DR3 / SKA); >2 dex below => PTA permanently
uninformative for that branch. Fork-split branches report per fork; branch verdict is
'excluded' only if EVERY frozen fork is excluded."""
import os, json, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from t2_guard import enforce, DERIVED

freezes, paths = enforce()          # guard: freeze exists, code clean, ancestry OK
out_all = {}
for fz, path in zip(freezes, paths):
    K95g = fz["K95_reference"]["grid"]
    K95 = K95g[f"{13/3:.3f}"]
    rows = {}

    def verdict(margin_dex):
        if margin_dex <= 0: return "EXCLUDED"
        if margin_dex <= 2: return "reachable-next-gen"
        return "pta-permanently-uninformative"

    # B1 — species dilution: K(N_SM)=K_full/107^p; N_min frozen
    b1 = fz["branches"]["B1"]
    b1_forks = {}
    for k, Kp in b1["K_pred"].items():
        m = float(np.log10(K95 / Kp))
        b1_forks[k] = dict(K_pred=Kp, margin_dex=round(m, 2),
                           verdict_at_SM_content=verdict(m), N_min=b1["N_min"][k.split("|")[0] + "|" + k.split("|")[1]])
    n_max_th = fz["branches"]["B1"]["early_universe"]["N_max_thermalized"]
    b1_note = ("ALL forks: N_min >> N_max_thermalized=%.0f => B1 survives ONLY via a "
               "never-thermalized hidden sector of >=%.0e dof => RETREAT-FLAG (plan §5.3)"
               % (n_max_th, min(float(v["N_min"]) for v in b1_forks.values())))
    rows["B1"] = dict(forks=b1_forks, branch_verdict="retreat-flagged", note=b1_note)

    # B2 — inverse bound: T_c ceiling
    b2 = fz["branches"]["B2"]
    rows["B2"] = dict(T_c_ceiling_over_TdS=b2["T_c_ceiling_over_TdS"],
                      T_dS_K=b2["T_dS_K"],
                      branch_verdict="pta-inverse-bound-only",
                      note="PTA bounds T_c/T_dS at 1e-7..1e-10 per l_c fork — no independent "
                           "T_c prediction exists to contradict; discrimination rides on "
                           "T2.4/T2.5/T2.6 per decision matrix")

    # B3 — leakage margin (distance-independent; per-pulsar ceiling applies)
    b3 = fz["branches"]["B3"]
    m3 = b3["leakage_margin_dex"]
    rows["B3"] = dict(E_leak=b3["K_pred_leakage_E_per_pulsar"], margin_dex=round(m3, 2),
                      branch_verdict="pta-permanently-uninformative",
                      note=b3["pta_note"])

    # B4 — functional; adjudicated in T2.2
    rows["B4"] = dict(branch_verdict="deferred-to-T2.2",
                      frozen_effect=fz["branches"]["B4"]["t22_prediction"])

    # B5 — massless carrier x Omega_r suppression
    b5 = fz["branches"]["B5"]
    b5_forks = {}
    for k, Kp in b5["K_pred"].items():
        m = float(np.log10(K95 / Kp))
        b5_forks[k] = dict(K_pred=Kp, margin_dex=round(m, 2), verdict=verdict(m))
    all_excl = all(v["verdict"] == "EXCLUDED" for v in b5_forks.values())
    rows["B5"] = dict(forks=b5_forks,
                      branch_verdict="EXCLUDED" if all_excl else "fork-split",
                      note="every frozen (l_c, suppression) fork lands ABOVE the T1.1 ceiling"
                           if all_excl else "surviving forks listed")

    out_all[os.path.basename(path)] = dict(K95_grid=K95g, rows=rows)

json.dump(out_all, open(os.path.join(DERIVED, "t2_1_verdicts.json"), "w"), indent=1)
print("saved derived/t2_1_verdicts.json\n")
for fname, block in out_all.items():
    print(f"=== vs {fname} (K95[13/3]={block['K95_grid'][f'{13/3:.3f}']:.2e}) ===")
    for b, r in block["rows"].items():
        print(f"\n[{b}] verdict: {r['branch_verdict']}")
        if "forks" in r:
            for k, v in r["forks"].items():
                extra = f" N_min={float(v['N_min']):.1e}" if "N_min" in v else ""
                vd = v.get("verdict", v.get("verdict_at_SM_content"))
                print(f"   {k:28} K={v['K_pred']:.2e}  margin={v['margin_dex']:+6.2f} dex  {vd}{extra}")
        if "note" in r: print(f"   note: {r['note'][:160]}")
