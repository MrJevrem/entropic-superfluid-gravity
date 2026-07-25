"""T3 pre-registration guard (T3_PLAN §1; semantics identical to t2_guard).

Refuses to run when: (a) no T3 freeze file; (b) analysis/ dirty in git;
(c) freeze commit not an ancestor of HEAD; (d) any staged comparison-data file
predates the freeze. Returns ALL T3 freeze versions, newest last.
"""
import json, os, subprocess, sys, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DERIVED = os.path.join(ROOT, "derived")

def _git(*args):
    return subprocess.run(["git", "-C", ROOT] + list(args), capture_output=True, text=True)

def enforce(staged=(), allow_dirty_paths=()):
    freezes = sorted(glob.glob(os.path.join(DERIVED, "T3_locked_predictions*.json")),
                     key=os.path.getmtime)
    if not freezes:
        sys.exit("T3 GUARD: no freeze file (derived/T3_locked_predictions*.json). "
                 "Run analysis/t3_predictions_freeze.py first.")
    dirty = [l for l in _git("status", "--porcelain", "--", "analysis/").stdout.splitlines()
             if l.strip() and not any(p in l for p in allow_dirty_paths)]
    if dirty:
        sys.exit("T3 GUARD: analysis/ dirty in git — commit before comparing:\n" + "\n".join(dirty))
    latest = json.load(open(freezes[-1]))
    commit = latest.get("provenance", {}).get("git_commit", "")
    if commit:
        if _git("merge-base", "--is-ancestor", commit, "HEAD").returncode != 0:
            sys.exit(f"T3 GUARD: freeze commit {commit[:8]} is not an ancestor of HEAD.")
    t_freeze = min(os.path.getmtime(f) for f in freezes)
    late = [s for s in staged if os.path.exists(s) and os.path.getmtime(s) < t_freeze]
    if late:
        sys.exit("T3 GUARD: comparison data was staged BEFORE the freeze (peeking risk): "
                 + ", ".join(late))
    return [json.load(open(f)) for f in freezes], freezes

if __name__ == "__main__":
    fz, paths = enforce()
    print(f"T3 guard OK: {len(fz)} freeze version(s): {[os.path.basename(p) for p in paths]}")
