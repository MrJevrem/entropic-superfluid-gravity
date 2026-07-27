#!/usr/bin/env python
"""
T3.9 — per-system stellar photometry for the Sun+09 groups (Paper V sec.7 (ii)).

Replaces the tabulated-LK stellar budget with independent 2MASS XSC K-band
photometry per group: Simbad-resolved centers, XSC cone within r500, luminosity
distances at the group z (H0=70, Om=0.3), M_K,sun = 3.28, k-corr -2.1z.
Two measured budgets bracket the systematics: BCG-only (brightest XSC source)
and summed bright members (no statistical background subtraction — biases M*
high, D low; logged). The t3_3 break fit is rerun with ONLY the stellar term
swapped; everything else (statistic, window, classes, weights, fit) frozen.
"""
import os, sys, re, json, io, zipfile, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import curve_fit
from t3_guard import enforce, ROOT, DERIVED

TEX = os.path.join(ROOT, "data/T3/sun09/ms.tex")
CACHE = os.path.join(ROOT, "data/T3/mixed/xsc_photometry.json")

# ---------- (A) acquire photometry (cached) ----------
def simbad(name):
    try:
        url = "https://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident=" + urllib.parse.quote(name)
        with urllib.request.urlopen(url, timeout=60) as r: t = r.read().decode(errors="ignore")
        m = re.search(r"Coordinates\(ICRS[^)]*\):\s*([\d]+)\s+([\d]+)\s+([\d.]+)\s+([+-][\d]+)\s+([\d]+)\s+([\d.]+)", t)
        if not m: return None
        ra = (float(m.group(1)) + float(m.group(2))/60 + float(m.group(3))/3600)*15
        d = float(m.group(4)); dec = abs(d) + float(m.group(5))/60 + float(m.group(6))/3600
        return ra, np.copysign(dec, d)
    except Exception: return None

def xsc(ra, dec, rs_arcmin):
    try:
        url = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=VII/233/xsc&"
               f"-c={ra:.5f}{dec:+.5f}&-c.rm={rs_arcmin:.2f}&-out.max=500")
        with urllib.request.urlopen(url, timeout=90) as r: t = r.read().decode(errors="ignore")
        lines = [l for l in t.splitlines() if l and not l.startswith("#")]
        hdr = next((l for l in lines if "K.ext" in l), None)
        if not hdr: return []
        ki = hdr.split("\t").index("K.ext")
        ks = []
        for l in lines:
            c = l.split("\t")
            if len(c) > ki:
                try: ks.append(float(c[ki]))
                except ValueError: pass
        return sorted(ks)
    except Exception: return []

if not os.path.exists(CACHE):
    tex0 = open(TEX).read()
    norm0 = lambda s: re.sub(r"\(.*?\)|[\s~]", "", s).replace("\\", "")
    names = []
    for line in tex0.split(r"\caption{The group sample")[1].split(r"\end{tabular}")[0].splitlines():
        cells = line.split("&")
        if len(cells) >= 8 and re.search(r"\d", cells[1]): names.append(norm0(cells[0]))
    # r500 per group for the cone radius
    r500s = {}
    for line in tex0.split(r"\caption{Derived properties of groups (I")[1].split(r"\end{tabular}")[0].splitlines():
        cells = line.split("&")
        if len(cells) >= 9:
            n = norm0(cells[0]); m = re.search(r"([\d.]+)", cells[3])
            if m: r500s[n] = float(m.group(1))
    def grab(name):
        c = simbad(name) or simbad(name.replace("RXJ", "RX J").replace("ESO", "ESO "))
        if not c: return name, None
        return name, dict(ra=c[0], dec=c[1])
    with ThreadPoolExecutor(6) as ex:
        coords = dict(ex.map(grab, names))
    out = {}
    for n, c in coords.items():
        if not c: out[n] = dict(fail="no coords"); continue
        out[n] = c
    json.dump(out, open(CACHE, "w"), indent=1)
    print(f"resolved {sum(1 for v in out.values() if 'ra' in v)}/{len(out)} group centers -> cache (cones next run)")

freezes, _ = enforce(staged=[TEX, CACHE], allow_dirty_paths=("t3_9_stellar_photometry.py",))
FZ = freezes[-1]
A0 = FZ["statistic"]["g_RAR_primary"]["a0_m_s2"]
GLO, GHI = FZ["statistic"]["g_bar_window_m_s2"]
H0 = FZ["statistic"]["cosmology"]["H0"]; OM = FZ["statistic"]["cosmology"]["Om"]
G = 6.674e-11; MSUN = 1.98892e30; KPC = 3.0857e19; C_KMS = 2.998e5
rho_c0 = 3 * (H0*1e3/(KPC*1e3))**2 / (8*np.pi*G)
g_rar = lambda gb: gb / (1 - np.exp(-np.sqrt(gb / A0)))
clean = lambda c: re.sub(r"\\pm|[\$\*{}]|\^\{[^}]*\}|_\{[^}]*\}|\\tablenotemark\{[a-z]\}|~", " ",
                         re.sub(r"\^\{[^}]*\}|_\{[^}]*\}", " ", c))
num1 = lambda c: (lambda m: float(m.group(1)) if m else np.nan)(re.search(r"([\d.]+)", clean(c)))
norm = lambda s: re.sub(r"\(.*?\)|[\s~]", "", s).replace("\\", "")

def dl_m(z):
    zs = np.linspace(0, z, 200)
    dc = np.trapezoid(1/np.sqrt(OM*(1+zs)**3 + 1-OM), zs)*C_KMS/H0
    return (1+z)*dc*1e3*KPC   # Mpc -> m

def lk_of(K, z):
    MK = K - 5*np.log10(dl_m(z)/ (10*3.0857e16)) + 2.1*z
    return 10**(-0.4*(MK - 3.28))

cache = json.load(open(CACHE))
tex = open(TEX).read()
zLK = {}
for line in tex.split(r"\caption{The group sample")[1].split(r"\end{tabular}")[0].splitlines():
    cells = line.split("&")
    if len(cells) >= 8 and re.search(r"\d", cells[1]):
        name = norm(cells[0]); z, LK = num1(cells[1]), num1(cells[7])
        if name and not np.isnan(z): zLK[name] = (z, LK)

# per-group photometric budgets (query cones now if missing)
phot = {}
for n, (z, LK) in zLK.items():
    c = cache.get(n, {})
    if "ra" not in c: continue
    if "ks" not in c:
        # r500 in kpc from derived table (reparse quickly)
        pass
    phot[n] = c
need = [n for n in phot if "ks" not in phot[n]]
if need:
    r500s = {}
    for line in tex.split(r"\caption{Derived properties of groups (I")[1].split(r"\end{tabular}")[0].splitlines():
        cells = line.split("&")
        if len(cells) >= 9:
            nn = norm(cells[0]); m = re.search(r"([\d.]+)", cells[3])
            if m: r500s[nn] = float(m.group(1))
    def cone(n):
        z = zLK[n][0]; r500 = r500s.get(n, 500.0)
        theta = min(np.degrees(r500*KPC/ (dl_m(z)/(1+z)**2))*60, 30.0)
        ks = xsc(phot[n]["ra"], phot[n]["dec"], theta)
        bg = xsc(phot[n]["ra"] + 2.5, phot[n]["dec"], theta)     # offset control cone (empirical background)
        return n, ks, bg
    with ThreadPoolExecutor(6) as ex:
        for n, ks, bg in ex.map(cone, need):
            phot[n]["ks"] = ks; phot[n]["ks_bg"] = bg
    json.dump({**cache, **phot}, open(CACHE, "w"), indent=1)
    print(f"XSC cones fetched for {len(need)} groups")

# ---------- (B) rebuild rows09 under three stellar budgets ----------
def build(budget):
    rows = []
    for line in tex.split(r"\caption{Derived properties of groups (I")[1].split(r"\end{tabular}")[0].splitlines():
        cells = line.split("&")
        if len(cells) < 9 or not re.search(r"\d", "".join(cells[1:3])): continue
        name = norm(cells[0])
        if name not in zLK: continue
        z, LK = zLK[name]
        T500 = num1(cells[1]); r500 = num1(cells[3]); r2500 = num1(cells[4])
        M500tab = num1(cells[5]); fg500 = num1(cells[6]); fg2500 = num1(cells[7])
        if np.isnan(T500) or np.isnan(r2500) or np.isnan(fg2500): continue
        rho_cz = rho_c0 * (OM*(1+z)**3 + 1-OM)
        if budget == "s09":
            Mstar = (10**LK)*MSUN if not np.isnan(LK) else 0.0
        else:
            ks = phot.get(name, {}).get("ks", [])
            if not ks: Mstar = (10**LK)*MSUN if not np.isnan(LK) else 0.0   # fallback, logged
            else:
                lks = [lk_of(k, z) for k in ks if k < 14.0]
                if budget == "bcg": L = max(lks) if lks else 0.0
                elif budget == "sum": L = sum(lks)
                else:                                                        # sum_bg: background-subtracted
                    bg = [lk_of(k, z) for k in phot.get(name, {}).get("ks_bg", []) if k < 14.0]
                    L = max(sum(lks) - sum(bg), max(lks) if lks else 0.0)
                Mstar = L*MSUN                                              # M/L_K = 1 (frozen)
        pts = []
        r_m = r2500*KPC; M2500 = (4*np.pi/3)*2500*rho_cz*r_m**3
        gobs = G*M2500/r_m**2; gbar = G*(fg2500*M2500 + 0.7*Mstar)/r_m**2
        if GLO <= gbar <= GHI: pts.append(np.log10(gobs) - np.log10(g_rar(gbar)))
        if not np.isnan(fg500) and not np.isnan(r500):
            r_m = r500*KPC
            M500 = M500tab*1e13*MSUN if not np.isnan(M500tab) else (4*np.pi/3)*500*rho_cz*r_m**3
            gobs5 = G*M500/r_m**2; gbar5 = G*(fg500*M500 + Mstar)/r_m**2
            if GLO <= gbar5 <= GHI: pts.append(np.log10(gobs5) - np.log10(g_rar(gbar5)))
        if pts:
            rows.append(dict(name=name, sigma=400.4*np.sqrt(T500), D=float(np.mean(pts)), Mstar=Mstar/MSUN))
    return rows

# other classes (identical to t3_3)
base_pts = []
zf = zipfile.ZipFile(os.path.join(ROOT, "data/T2/sparc/Rotmod_LTG.zip"))
for nm in zf.namelist():
    if not nm.endswith("_rotmod.dat"): continue
    try: arr = np.loadtxt(io.BytesIO(zf.read(nm)), comments="#")
    except Exception: continue
    if arr.ndim != 2 or arr.shape[0] < 6 or arr.shape[1] < 6: continue
    r, vobs, _, vgas, vdisk, vbul = (arr[:, i] for i in range(6))
    r_m = r*KPC; km = 1e3
    gobs = (vobs*km)**2/r_m
    gbar = (vgas*km*np.abs(vgas*km) + 0.5*(vdisk*km)**2 + 0.7*(vbul*km)**2)/r_m
    m = (gbar > 0) & (gbar >= GLO) & (gbar <= GHI) & (gobs > 0)
    if m.sum() >= 3:
        base_pts.append((float(np.mean(vobs[-3:]))/np.sqrt(2),
                         float(np.mean(np.log10(gobs[m]) - np.log10(g_rar(gbar[m])))), "SPARC"))
for r_ in json.load(open(os.path.join(DERIVED, "t3_2_shelf.json")))["rows"]:
    if r_["rho_env"] < 0.5: base_pts.append((r_["sigma"], r_["D"], "ETG_field"))
v06 = [(r_["sigma"], r_["D"], "V06") for r_ in json.load(open(os.path.join(DERIVED, "t3_1_anchors.json")))["v06"]]

CLS_SIG = {"SPARC": .20, "ETG_field": .22, "Sun09": .12, "V06": .08}
def model(lnS, a, b, lnSb, wid): return a + b/(1 + np.exp(-(lnS - lnSb)/wid))
def fit_budget(budget):
    rows = build(budget)
    pts = base_pts + [(r["sigma"], r["D"], "Sun09") for r in rows] + v06
    sig = np.array([p[0] for p in pts]); D = np.array([p[1] for p in pts])
    w = np.array([CLS_SIG[p[2]] for p in pts])
    p, _ = curve_fit(model, np.log(sig), D, p0=[0., .45, np.log(500.), .3],
                     sigma=w, bounds=([-.2, 0, np.log(100), .02], [.2, 1, np.log(2000), 1.5]), maxfev=20000)
    s09D = np.array([r["D"] for r in rows])
    return rows, float(np.exp(p[2])), float(s09D.mean())

print("=== T3.9: the break fit under measured per-system photometry ===")
res = {}
for budget, lab in (("s09", "S09 tabulated LK (baseline)"), ("bcg", "2MASS XSC BCG-only"), ("sum", "2MASS XSC summed (raw)"), ("sum_bg", "2MASS XSC summed, bg-subtracted")):
    rows, sb, dmean = fit_budget(budget)
    nph = sum(1 for r in rows if phot.get(r["name"], {}).get("ks"))
    res[budget] = dict(sigma_b=sb, S09_meanD=dmean, n=len(rows), n_phot=nph)
    print(f"{lab:32s}: sigma_b = {sb:5.0f} km/s   <D_S09> = {dmean:+.3f}   (N={len(rows)}, photometric {nph})")

# photometric vs tabulated light comparison
comp = []
for n, (z, LK) in zLK.items():
    ks = phot.get(n, {}).get("ks", [])
    if ks and not np.isnan(LK):
        lsum = sum(lk_of(k, z) for k in ks if k < 14.0)
        if lsum > 0: comp.append(np.log10(lsum) - LK)
comp = np.array(comp)
print(f"\nphotometric(sum) vs S09 tabulated log LK: offset {comp.mean():+.2f} dex, scatter {comp.std():.2f} dex (N={len(comp)})")
print(f"stellar-budget fork before: sigma_b extremes [215, 352] (conventions); now: [{min(r['sigma_b'] for r in res.values()):.0f}, {max(r['sigma_b'] for r in res.values()):.0f}] (measured band)")

json.dump(dict(frozen_against=FZ["provenance"], budgets=res, phot_vs_tab=dict(offset=float(comp.mean()), scatter=float(comp.std()), n=len(comp)),
               conventions=["XSC Kmag < 14, cone r500 (<=30'), M_Ksun=3.28, kcorr -2.1z, M/L_K=1 frozen",
                            "no statistical background subtraction (biases summed M* high, D low; logged)",
                            "groups without photometry fall back to S09 LK (logged per count)"]),
          open(os.path.join(DERIVED, "t3_9_stellar_photometry.json"), "w"), indent=1)
print("wrote derived/t3_9_stellar_photometry.json")
