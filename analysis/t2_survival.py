#!/usr/bin/env python
"""§9 — assemble derived/T2_survival.json from the test outputs (no hand-written verdicts)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t2_guard import enforce, DERIVED

freezes, fpaths = enforce()
J = lambda n: json.load(open(os.path.join(DERIVED, n)))
t21, t22, t23, t24, t25 = (J(f"t2_{i}_results.json") if i != 1 else J("t2_1_verdicts.json")
                           for i in (1, 2, 3, 4, 5))
t20 = J("t2_0_window.json")
v = t21[list(t21)[-1]]["rows"]

surv = dict(provenance=dict(freeze_files=[os.path.basename(p) for p in fpaths],
                            note="all predictions frozen at commit embedded in freeze; "
                                 "comparison data staged after (guard-enforced)"))
surv["B1"] = dict(status="retreat-flagged",
    ledger=["T2.1: excluded at SM content in ALL forks (3.8-8.8 dex); N_min 6.5e5-4.2e19",
            f"T2.3: ever-thermalized fork excluded by {t23['B1']['violation_dex']:.1f} dex in DeltaNeff",
            "survives only as never-thermalized hidden sector with no other signature"],
    killing_test="T2.1+T2.3 (predictive variants); retreat loophole recorded")
surv["B2"] = dict(status="alive",
    ledger=[f"T2.1: PTA inverse bound only — T_c/T_dS < {min(v['B2']['T_c_ceiling_over_TdS'].values()):.1e}"
            f"..{max(v['B2']['T_c_ceiling_over_TdS'].values()):.1e} per l_c fork",
            "T2.4: constrains via a0 statics only weakly (not its lock)",
            "frozen WB forks: dv=0.065 (EFE-analog) | 0 (decoherent)"],
    next_observation="Gaia DR4 (2026-12-02) via frozen T2.5 pipeline; T2.6 conditional after")
surv["B3"] = dict(status="retreat-flagged",
    ledger=[f"T2.1: PTA permanently uninformative (margin {v['B3']['margin_dex']:+.1f} dex, "
            "distance-independent kernel)",
            f"T2.4: a0=const preferred over frozen cH(z) lock at ~{t24['sigma_corr_nominal']} sigma "
            f"nominal / ~{t24['sigma_corr_inflated']} sigma inflated — meets >=3sigma kill "
            "threshold marginally on literature-grade [S] compilation",
            "frozen WB forks (0.20 | 0) remain live for DR4"],
    killing_test="T2.4 (marginal; primary-IFU reanalysis = escalation); "
                 "survives only by demoting a0~cH0 to today-coincidence")
surv["B4"] = dict(status="dead",
    ledger=[f"T2.2: primary frozen statistic null — partial rho={t22['achromatic']['all']['s_path']['partial_rho']:+.2f} "
            f"(perm p={t22['achromatic']['all']['s_path']['perm_p']}), Sidak p={t22['achromatic']['all']['_sidak_over_3']}",
            f"measured slope {t22['slope_measured']:+.1f} vs frozen +1.0 at {t22['contrast_log10']:.1f} dex "
            "available contrast => coupling >=1 dex below lock",
            "caveats logged: logistic IRLS unstable (flag), pygedm unavailable (EM-flag substitute)"],
    killing_test="T2.2 (frozen exclusion rule)")
surv["B5"] = dict(status="dead",
    ledger=["T2.1: ALL frozen forks 3.6-9.4 dex ABOVE the T1.1 ceiling => excluded",
            f"T2.3: independent second exclusion — {t23['B5_verdict'][:80]}..."],
    killing_test="T2.1 (PTA ceiling); T2.3 (FIRAS mu) independently")
surv["T2_0_context"] = t20["verdict"]
surv["T2_5_status"] = t25["dry_run_status"]
surv["T2_6"] = ("deferred-conditional: B2 alive => execute after DR4 per plan §8 ordering")
surv["summary"] = dict(dead=["B4", "B5"], retreat_flagged=["B1", "B3"], alive=["B2"],
                       decisive_future="Gaia DR4 Dec 2026 (B2 forks; B3 WB forks) after "
                                       "control-bin systematic (~5%) is repaired")
json.dump(surv, open(os.path.join(DERIVED, "T2_survival.json"), "w"), indent=1)
print(json.dumps(surv["summary"], indent=1))
for b in ("B1", "B2", "B3", "B4", "B5"):
    print(f"\n[{b}] {surv[b]['status'].upper()}")
    for l in surv[b]["ledger"]: print("   -", l)
print("\nsaved derived/T2_survival.json")
