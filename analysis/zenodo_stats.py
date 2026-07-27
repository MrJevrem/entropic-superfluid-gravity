#!/usr/bin/env python
"""
Snapshot Zenodo + GitHub usage stats with a UTC timestamp, appending to
derived/usage_stats_log.jsonl. Re-run any time (cron-able); view INCREMENTS
between snapshots can then be correlated with outreach events (emails, posts) —
the closest ethical proxy for "did anyone look" that the platforms allow.
"""
import json, subprocess, urllib.request, datetime, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t3_guard import DERIVED

RECS = {"ms_v1": 21561507, "ms_v2": 21571242, "ms_v3": 21572393, "ms_v4": 21626776, "ms_v5": 21627204, "ms_v6": 21627781,
        "code_v1.0": 21560681, "code_v1.1": 21571552, "code_v1.2": 21626861}
snap = {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
for name, rid in RECS.items():
    try:
        with urllib.request.urlopen(f"https://zenodo.org/api/records/{rid}", timeout=60) as r:
            s = json.loads(r.read().decode(), strict=False).get("stats", {})
        snap[name] = {k: s.get(k, 0) for k in ("views", "unique_views", "downloads", "unique_downloads",
                                                "version_views", "version_unique_views",
                                                "version_downloads", "version_unique_downloads")}
        # NB Zenodo semantics: plain keys = ALL-VERSIONS concept aggregate (identical across
        # sibling versions); version_* keys = THIS record's own page. REST API hits are
        # filtered by their COUNTER pipeline and do not count as views (verified 2026-07-27).
    except Exception as e:
        snap[name] = {"error": str(e)[:60]}
try:
    out = subprocess.run(["gh", "api", "repos/MrJevrem/entropic-superfluid-gravity/traffic/popular/referrers",
                          "--jq", "[.[] | {r: .referrer, c: .count, u: .uniques}]"], capture_output=True, text=True, timeout=60)
    snap["gh_referrers"] = json.loads(out.stdout) if out.returncode == 0 and out.stdout.strip() else []
except Exception:
    snap["gh_referrers"] = []
for ep, key in (("views", "gh_views"), ("clones", "gh_clones")):
    try:
        out = subprocess.run(["gh", "api", f"repos/MrJevrem/entropic-superfluid-gravity/traffic/{ep}",
                              "--jq", "{count: .count, uniques: .uniques}"], capture_output=True, text=True, timeout=60)
        snap[key] = json.loads(out.stdout) if out.returncode == 0 else {"error": out.stderr[:60]}
    except Exception as e:
        snap[key] = {"error": str(e)[:60]}

log = os.path.join(DERIVED, "usage_stats_log.jsonl")
with open(log, "a") as f:
    f.write(json.dumps(snap) + "\n")
prev = None
lines = open(log).read().strip().splitlines()
if len(lines) >= 2:
    prev = json.loads(lines[-2])
print(f"snapshot @ {snap['utc']}")
for name in RECS:
    v = snap[name].get("unique_views", "?"); vv = snap[name].get("version_unique_views", "?")
    d = f" (+{snap[name]['unique_views'] - prev[name]['unique_views']})" if prev and "unique_views" in prev.get(name, {}) else ""
    print(f"  {name:10s} concept uniq: {v}{d}  own-page uniq: {vv}  downloads: {snap[name].get('unique_downloads', '?')}")
print(f"  github: views {snap.get('gh_views', {})}, clones {snap.get('gh_clones', {})}")
print(f"appended to {log} ({len(lines)} snapshots)")
