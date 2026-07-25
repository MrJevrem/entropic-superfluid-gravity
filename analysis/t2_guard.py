"""T2 pre-registration guard (T2_PLAN §1). Comparison scripts import and call `enforce()`.

Refuses to run when:
  (a) the freeze file is missing;
  (b) analysis code is dirty in git (uncommitted changes to analysis/*.py);
  (c) the freeze file's embedded commit is not an ancestor of HEAD;
  (d) any staged comparison-data file passed in `staged` predates the freeze — i.e. the
      freeze must have been created BEFORE the comparison data was staged (no peeking).
Multiple freeze versions: returns ALL freeze files (T2_locked_predictions*.json), newest last;
comparisons must report against every version (§1.3).
"""
import json, os, subprocess, sys, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DERIVED = os.path.join(ROOT, "derived")

def _git(*args):
    return subprocess.run(["git", "-C", ROOT] + list(args), capture_output=True, text=True)

def enforce(staged=(), allow_dirty_paths=()):
    freezes = sorted(glob.glob(os.path.join(DERIVED, "T2_locked_predictions*.json")),
                     key=os.path.getmtime)
    if not freezes:
        sys.exit("T2 GUARD: no freeze file (derived/T2_locked_predictions*.json). "
                 "Run analysis/t2_predictions_freeze.py first.")
    dirty = [l for l in _git("status", "--porcelain", "--", "analysis/").stdout.splitlines()
             if l.strip() and not any(p in l for p in allow_dirty_paths)]
    if dirty:
        sys.exit(f"T2 GUARD: analysis/ dirty in git — commit before comparing:\n" + "\n".join(dirty))
    latest = json.load(open(freezes[-1]))
    commit = latest.get("provenance", {}).get("git_commit", "")
    if commit:
        ok = _git("merge-base", "--is-ancestor", commit, "HEAD").returncode == 0
        if not ok:
            sys.exit(f"T2 GUARD: freeze commit {commit[:8]} is not an ancestor of HEAD.")
    t_freeze = min(os.path.getmtime(f) for f in freezes)
    late = [s for s in staged if os.path.exists(s) and os.path.getmtime(s) < t_freeze]
    # NOTE the inequality direction: staged data must be YOUNGER than the freeze.
    if late:
        sys.exit("T2 GUARD: comparison data was staged BEFORE the freeze (peeking risk): "
                 + ", ".join(late))
    return [json.load(open(f)) for f in freezes], freezes

if __name__ == "__main__":
    fz, paths = enforce()
    print(f"guard OK: {len(fz)} freeze version(s): {[os.path.basename(p) for p in paths]}")
