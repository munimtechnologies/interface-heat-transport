# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_cover 0.1.0: the conserving "cover" band rule for the no-fit multiband two-stream model
(PLAN amendment 0.5.11 section 4, tests V0-V5; the case the rule did not foresee is specified in
amendment 0.5.12 section 2, written before this module ran).

tbc-fill-v1 showed that the 5-6x deficit of TiN interfaces below their own DMM came mostly from
fine frequency bins (0.2 THz) that have a radiation conductance K > 0 on the dense harmonic mesh
but no mode on the 16^3 lifetime mesh: tbc_bte.lengths_from_kappa gives them interpolated
lengths, so the model's TiN conductivity was 132 instead of 44.6 W/m/K. Merging the bins below
the first jointly sampled bin fixed most of it, but TiN's filled bins reach 3.2 THz and stay
filled where merging stopped early (Si, Ge).

Cover rule. From 0 THz, consecutive fine bins form one band until, on EACH side, the band holds a
bin sampled on that side's lifetime mesh (K > 0 and kappa > 0) or holds no state at all on that
side (K = 0 in every bin: nothing to fill; amendment 0.5.12 sec. 2). A side's bins above its
last sampled bin join the band that holds that last sampled bin (the pre-registered "tail joins
the last band", applied per side so that one side's spectrum above the other's cutoff is not
folded onto it). In each band K and kappa are summed per side; the length is
kappa_band / (2 K_band); the DMM is recomputed on the same bands (tbc_bte does that from the
banded K). No length is interpolated, so each side's model conductivity equals its binned
phono3py value by construction.

Also reported (never graded): the LITERAL reading of 0.5.11 sec. 4 (close a band only when both
sides hold a sampled bin; any tail joins the last band), base and merge (tbc-fill-v1 rules,
closure), and the four earlier pairs at 300 K. tbc_bte, tbc_anh and tbc_fill are imported
unchanged. Last output line `COV:`.
"""
VERSION = "mv2 tbc-cover 0.1.0"

PAIRS = (("TiN", "C"), ("TiN", "Si"), ("TiN", "Ge"))
T_HIGH = {"TiN/C": 900, "TiN/Si": 700, "TiN/Ge": 700}
NODES = 400
MODEL = "orb-v3"

# tbc-fill-v1 (results/tin/fill/fill-summary.json): G base, G merge (MW m^-2 K^-1) and the upper
# edge of the merged band (THz)
FILL = {"TiN/C": {300: (43.251743, 192.927257, 4.8), 900: (52.804850, 235.766812, 4.8)},
        "TiN/Si": {300: (52.342830, 316.966128, 1.0), 700: (61.473291, 357.393892, 1.0)},
        "TiN/Ge": {300: (51.474463, 169.585272, 0.6), 700: (59.337807, 187.995254, 0.6)}}
G_INEL = {"TiN/C": {300: 5.7, 900: 44.1}, "TiN/Si": {300: 56.4, 700: 172.4},
          "TiN/Ge": {300: 61.8, 700: 188.2}}
MEASURED = {"TiN/C": {300: 449.0, 900: 910.0}, "TiN/Si": {300: 507.0, 700: 916.0},
            "TiN/Ge": {300: 190.0, 700: 281.0}}

# earlier pairs (300 K, orb-v3 npz of the stage that solved them; oriented geometry as there:
# M_or with the c-axis conductivity for wurtzite)
EARLIER = (("Al/Si", "results/anh/v1", "Al", "Si"),
           ("Al/GaN", "results/anh/v1", "Al", "GaN"),
           ("AlN/GaN", "results/anh/v1", "AlN", "GaN"),
           ("ZnO/GaN", "results/nonmetal/zno-v1", "ZnO", "GaN"))
KDIR = {"AlN": "z", "GaN": "z", "ZnO": "z"}

# pre-registered thresholds (amendment 0.5.11 section 4)
V0_ABS, V0_REL = 0.5, 0.01
V1_REL = 0.01
V2_C_REL = 0.25
V3_LO, V3_HI = 0.4, 1.4
V4_REL = 0.10
V4_MIN_BINS = 2
V5_FACTOR = 1.2


# ----------------------------------------------------------------------------- bands
def _sampled(K, kap, k):
    return K[k] > 0.0 and kap[k] > 0.0


def last_sampled(K, kap):
    idx = [k for k in range(len(K)) if _sampled(K, kap, k)]
    return idx[-1] if idx else None


def cover_bands(Ka, kapa, Kb, kapb, min_bins=1, literal=False):
    """Bands as inclusive (start, end) index pairs covering every fine bin.

    literal=False: the rule of the module docstring (empty side counts as covered; a side's
    bins above its last sampled bin are ignored for closing and joined later, see banded()).
    literal=True: close only when both sides hold a sampled bin; the tail joins the last band."""
    n = len(Ka)
    sides = ((Ka, kapa), (Kb, kapb))
    lasts = [last_sampled(K, kap) for K, kap in sides]
    if lasts[0] is None or lasts[1] is None:
        raise ValueError("a side has no sampled bin")
    bands = []
    start = 0
    for i in range(n):
        if i - start + 1 < min_bins:
            continue
        ok = True
        for (K, kap), ls in zip(sides, lasts):
            has = any(_sampled(K, kap, k) for k in range(start, i + 1))
            if literal:
                ok = ok and has
            else:
                empty = all(not (K[k] > 0.0) or k > ls for k in range(start, i + 1))
                ok = ok and (has or empty)
        if ok:
            bands.append((start, i))
            start = i + 1
    if start < n:
        if bands:
            bands[-1] = (bands[-1][0], n - 1)
        else:
            bands.append((0, n - 1))
    return bands


def banded(K, kap, bands, per_side_tail=True):
    """Sum K and kappa into bands. With per_side_tail, this side's bins above its last sampled
    bin are added to the band that holds that bin. Returns (K_band, kappa_band) lists."""
    K = [float(x) for x in K]
    kap = [float(x) for x in kap]
    ls = last_sampled(K, kap)
    Kb, kb = [], []
    home = None
    for j, (s, e) in enumerate(bands):
        hi = e if not per_side_tail else min(e, ls)
        Kb.append(sum(K[s:hi + 1]) if hi >= s else 0.0)
        kb.append(sum(kap[s:hi + 1]) if hi >= s else 0.0)
        if s <= ls <= e:
            home = j
    if per_side_tail and ls < len(K) - 1:
        Kb[home] += sum(K[ls + 1:])
        kb[home] += sum(kap[ls + 1:])
    return Kb, kb


def band_widths(bands, edges_thz):
    return [edges_thz[e + 1] - edges_thz[s] for s, e in bands]


def within(v, ref, abs_tol, rel_tol):
    return abs(v - ref) <= max(abs_tol, rel_tol * abs(ref))


# ----------------------------------------------------------------------------- grades
def grade(rows, tests_ok=True, earlier=None):
    """rows[pair][T] dicts with base, merge, first_band_hi_thz, cover, dmm_cover, kap1_cover,
    kap2_cover, kap1_bins, kap2_bins, cover_coarse. earlier: {pair: {base, cover}}."""
    g = {}
    bad = []
    for p, by in FILL.items():
        for T, (gb, gm, hi) in by.items():
            r = rows.get(p, {}).get(T)
            if r is None:
                bad.append(f"{p} {T} missing")
                continue
            if abs(r["first_band_hi_thz"] - hi) > 1e-6:
                bad.append(f"{p} {T} first band to {r['first_band_hi_thz']:.1f} vs {hi:.1f} THz")
            if not within(r["base"], gb, V0_ABS, V0_REL):
                bad.append(f"{p} {T} base {r['base']:.1f} vs {gb:.1f}")
            if not within(r["merge"], gm, V0_ABS, V0_REL):
                bad.append(f"{p} {T} merge {r['merge']:.1f} vs {gm:.1f}")
    if not tests_ok:
        bad.append("unit tests")
    g["V0"] = "PASS" if not bad else "FAIL(" + "; ".join(bad) + ")"
    allr = [(p, T, r) for p, by in rows.items() for T, r in by.items()]
    if not allr:
        for k in ("V1", "V2", "V3", "V4", "V5"):
            g[k] = "UNDETERMINED(no rows)"
        return g
    b1 = []
    for p, T, r in allr:
        for s in ("1", "2"):
            m, ref = r["kap" + s + "_cover"], r["kap" + s + "_bins"]
            if not abs(m / ref - 1.0) <= V1_REL:
                b1.append(f"{p} {T} side{s} {m:.2f} vs {ref:.2f}")
    g["V1"] = "PASS" if not b1 else "FAIL(" + ", ".join(b1) + ")"
    b2 = []
    for p, T, r in allr:
        q = r["cover"] / r["merge"]
        if p == "TiN/C":
            if not abs(q - 1.0) <= V2_C_REL:
                b2.append(f"{p} {T} cover/merge {q:.3f}")
        elif not r["cover"] >= r["merge"]:
            b2.append(f"{p} {T} cover/merge {q:.3f}")
    g["V2"] = ("PASS" if not b2 else "FAIL(" + ", ".join(b2) + ")") + "(" + ", ".join(
        f"{p} {T} {r['cover'] / r['merge']:.3f}" for p, T, r in allr) + ")"
    b3 = [f"{p} {T} {r['cover'] / r['dmm_cover']:.3f}" for p, T, r in allr
          if not V3_LO <= r["cover"] / r["dmm_cover"] <= V3_HI]
    g["V3"] = ("PASS" if not b3 else "FAIL(" + ", ".join(b3) + ")") + "(" + ", ".join(
        f"{p} {T} {r['cover'] / r['dmm_cover']:.3f}" for p, T, r in allr) + ")"
    b4 = [f"{p} {T} {r['cover_coarse'] / r['cover']:.3f}" for p, T, r in allr
          if not abs(r["cover_coarse"] / r["cover"] - 1.0) <= V4_REL]
    g["V4"] = ("PASS" if not b4 else "FAIL(" + ", ".join(b4) + ")") + "(" + ", ".join(
        f"{p} {T} {r['cover_coarse'] / r['cover']:.3f}" for p, T, r in allr) + ")"
    z = (earlier or {}).get("ZnO/GaN")
    if not z or not z.get("base") or z.get("cover") is None:
        g["V5"] = "UNDETERMINED(ZnO/GaN not solved)"
    else:
        q = z["cover"] / z["base"]
        g["V5"] = ("PASS" if q >= V5_FACTOR else "FAIL") + f"(ZnO/GaN cover/base {q:.3f})"
    return g


# ----------------------------------------------------------------------------- driver
def main():
    import argparse
    import json
    import os
    import time

    import numpy as np

    import tbc_anh as ta
    import tbc_bte as tb
    import tbc_fill as tf

    ap = argparse.ArgumentParser()
    ap.add_argument("--npz-dir", required=True, help="root with results/... npz paths")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tests-ok", default="1")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    lines = []

    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"{VERSION} | {ta.VERSION} | {tb.VERSION} | {tf.VERSION}")

    def load(stage, m):
        z = np.load(os.path.join(args.npz_dir, stage, f"anh-{MODEL}-{m}.npz"))
        r = {k: z[k] for k in z.files}
        temps = [float(x) for x in r["temps"]]
        it = temps.index(300.0)
        r1 = dict(r)
        r1["vcell"] = 1.0
        s = float(ta.kappa_bins(r1, it, 300.0, r["edges_c"], "iso", with_isotope=True).sum())
        k_ph3 = float(np.mean(r["kappa_ph3"][it][:3]))
        r["vcell"] = s / k_ph3
        return r, temps

    def bins(mat, T, dirn="iso"):
        r, temps = mat
        i = temps.index(float(T))
        K = ta.K_bins(r["M_or"], r["edges_f"], r["edges_c"], float(T), classical=False)
        kap = ta.kappa_bins(r, i, float(T), r["edges_c"], dirn, with_isotope=True)
        return np.asarray(K, float), np.asarray(kap, float)

    def solve_cover(Ka, kapa, Kb, kapb, min_bins=1, literal=False):
        bands = cover_bands(list(Ka), list(kapa), list(Kb), list(kapb), min_bins=min_bins,
                            literal=literal)
        pst = not literal
        Ka2, ka2 = banded(Ka, kapa, bands, per_side_tail=pst)
        Kb2, kb2 = banded(Kb, kapb, bands, per_side_tail=pst)
        Ka2, ka2, Kb2, kb2 = (np.array(x) for x in (Ka2, ka2, Kb2, kb2))
        la, fa = tb.lengths_from_kappa(ka2, Ka2)
        lb, fb = tb.lengths_from_kappa(kb2, Kb2)
        s = tb.solve_two_stream(Ka2, Kb2, la, lb, n_nodes=NODES)
        s["bands"] = bands
        s["filled"] = [fa, fb]
        # kappa in bands where that side has no K: lost to the model (reported)
        s["kappa_lost"] = [float(np.sum(ka2[Ka2 <= 0.0])), float(np.sum(kb2[Kb2 <= 0.0]))]
        return s

    def solve_all(A, B, Ka, kapa, Kb, kapb, edges_thz):
        la, _ = tb.lengths_from_kappa(kapa, Ka)
        lb, _ = tb.lengths_from_kappa(kapb, Kb)
        base = tb.solve_two_stream(Ka, Kb, la, lb, n_nodes=NODES)
        Ka3, kapa3, Kb3, kapb3, j = tf.merge_low(list(Ka), list(kapa), list(Kb), list(kapb))
        Ka3, kapa3, Kb3, kapb3 = (np.array(x) for x in (Ka3, kapa3, Kb3, kapb3))
        la3, _ = tb.lengths_from_kappa(kapa3, Ka3)
        lb3, _ = tb.lengths_from_kappa(kapb3, Kb3)
        mg = tb.solve_two_stream(Ka3, Kb3, la3, lb3, n_nodes=NODES)
        cv = solve_cover(Ka, kapa, Kb, kapb)
        cc = solve_cover(Ka, kapa, Kb, kapb, min_bins=V4_MIN_BINS)
        lit = solve_cover(Ka, kapa, Kb, kapb, literal=True)
        w = band_widths(cv["bands"], edges_thz)
        r = {"base": base["G"] / 1e6, "dmm_base": base["G_dmm"] / 1e6, "merge": mg["G"] / 1e6,
             "dmm_merge": mg["G_dmm"] / 1e6, "merge_upto_thz": float(edges_thz[j + 1]),
             "cover": cv["G"] / 1e6, "dmm_cover": cv["G_dmm"] / 1e6,
             "simons_cover": cv["G_simons"] / 1e6, "linbath_cover": cv["G_linear_bath"] / 1e6,
             "cover_coarse": cc["G"] / 1e6, "dmm_coarse": cc["G_dmm"] / 1e6,
             "literal": lit["G"] / 1e6, "dmm_literal": lit["G_dmm"] / 1e6,
             "kap1_cover": cv["kappa1"], "kap2_cover": cv["kappa2"],
             "kap1_bins": float(kapa.sum()), "kap2_bins": float(kapb.sum()),
             "kap1_base": base["kappa1"], "kap2_base": base["kappa2"],
             "kappa_lost": cv["kappa_lost"], "filled_cover": cv["filled"],
             "capped_share": [cv["capped_share1"], cv["capped_share2"]],
             "n_bands": len(cv["bands"]), "n_bands_coarse": len(cc["bands"]),
             "n_bands_literal": len(lit["bands"]),
             "first_band_hi_thz": float(edges_thz[cv["bands"][0][1] + 1]),
             "band_width_min_max_thz": [min(w), max(w)],
             "wide_bands_thz": [[round(float(edges_thz[s]), 2), round(float(edges_thz[e + 1]), 2)]
                                for (s, e), x in zip(cv["bands"], w) if x > 0.2 + 1e-9],
             "last_sampled_thz": [float(edges_thz[last_sampled(list(Ka), list(kapa)) + 1]),
                                  float(edges_thz[last_sampled(list(Kb), list(kapb)) + 1])],
             "residual": max(x["residual"] for x in (base, mg, cv, cc, lit))}
        return r

    rows = {}
    mats = {m: load("results/tin/orb-v3", m) for m in ("TiN", "C", "Si", "Ge")}
    edges = np.asarray(mats["TiN"][0]["edges_c"], dtype=float)
    for m, mat in mats.items():
        if not np.allclose(edges, mat[0]["edges_c"]):
            raise SystemExit(f"edges_c differ for {m}")
    e_thz = edges / 1e12
    log(f"fine bin width {e_thz[1] - e_thz[0]:.3f} THz, {len(edges) - 1} bins to {e_thz[-1]:.1f} THz")
    for A, B in PAIRS:
        pair = f"{A}/{B}"
        rows[pair] = {}
        for T in (300, T_HIGH[pair]):
            Ka, kapa = bins(mats[A], T)
            Kb, kapb = bins(mats[B], T)
            r = solve_all(A, B, Ka, kapa, Kb, kapb, e_thz)
            r["inel"] = G_INEL[pair][T]
            r["measured"] = MEASURED[pair][T]
            r["total_cover"] = r["cover"] + r["inel"]
            rows[pair][T] = r
            log(f"COVER {pair} {T} K | first band 0-{r['first_band_hi_thz']:.1f} THz (merge "
                f"0-{r['merge_upto_thz']:.1f}) | G base {r['base']:.1f} merge {r['merge']:.1f} "
                f"cover {r['cover']:.1f} coarse {r['cover_coarse']:.1f} literal {r['literal']:.1f}"
                f" | DMM base {r['dmm_base']:.1f} merge {r['dmm_merge']:.1f} cover "
                f"{r['dmm_cover']:.1f} coarse {r['dmm_coarse']:.1f} literal {r['dmm_literal']:.1f}"
                f" | simons_c {r['simons_cover']:.1f} linbath_c {r['linbath_cover']:.1f} | kappa "
                f"cover/bins {r['kap1_cover']:.2f}/{r['kap1_bins']:.2f} {r['kap2_cover']:.2f}/"
                f"{r['kap2_bins']:.2f} lost {r['kappa_lost'][0]:.3f}/{r['kappa_lost'][1]:.3f} "
                f"filled {r['filled_cover']} capped {r['capped_share'][0]:.3f}/"
                f"{r['capped_share'][1]:.3f} | bands {r['n_bands']} (coarse "
                f"{r['n_bands_coarse']}, literal {r['n_bands_literal']}) width "
                f"{r['band_width_min_max_thz'][0]:.1f}-{r['band_width_min_max_thz'][1]:.1f} THz | "
                f"last sampled {r['last_sampled_thz'][0]:.1f}/{r['last_sampled_thz'][1]:.1f} THz"
                f" | info cover+inel {r['total_cover']:.0f} vs measured {r['measured']:.0f} | "
                f"resid {r['residual']:.1e}")
            log(f"WIDE {pair} {T} K bands wider than one bin (THz): {r['wide_bands_thz']}")
    del mats

    earlier = {}
    for pair, stage, A, B in EARLIER:
        try:
            ma, mb = load(stage, A), load(stage, B)
            ea = np.asarray(ma[0]["edges_c"], float)
            if not np.allclose(ea, mb[0]["edges_c"]):
                raise ValueError("edges_c differ")
            Ka, kapa = bins(ma, 300, KDIR.get(A, "iso"))
            Kb, kapb = bins(mb, 300, KDIR.get(B, "iso"))
            r = solve_all(A, B, Ka, kapa, Kb, kapb, ea / 1e12)
            earlier[pair] = r
            log(f"EARLIER {pair} 300 K ({stage}, {MODEL}, oriented) | G base {r['base']:.1f} "
                f"merge {r['merge']:.1f} cover {r['cover']:.1f} coarse {r['cover_coarse']:.1f} "
                f"literal {r['literal']:.1f} | DMM base {r['dmm_base']:.1f} cover "
                f"{r['dmm_cover']:.1f} | kappa cover/bins {r['kap1_cover']:.2f}/"
                f"{r['kap1_bins']:.2f} {r['kap2_cover']:.2f}/{r['kap2_bins']:.2f} | bands "
                f"{r['n_bands']} | first band 0-{r['first_band_hi_thz']:.1f} THz")
            del ma, mb
        except Exception as e:  # noqa: BLE001
            earlier[pair] = {"error": f"{type(e).__name__}: {e}"}
            log(f"EARLIER {pair}: FAILED {type(e).__name__}: {e}")

    g = grade(rows, tests_ok=args.tests_ok == "1", earlier=earlier)
    for k, v in g.items():
        log(f"{k}: {v}")
    with open(os.path.join(args.out, "cover-summary.json"), "w") as fh:
        json.dump({"version": VERSION, "grades": g,
                   "rows": {p: {str(T): v for T, v in by.items()} for p, by in rows.items()},
                   "earlier": earlier, "seconds": round(time.time() - t0, 1)}, fh, indent=1,
                  default=float)
    log("COV: " + " | ".join(f"{k} {str(v).split('(')[0]}" for k, v in g.items()) + " | " +
        " ; ".join(f"{p} {T} base/merge/cover/dmm_c {r['base']:.0f}/{r['merge']:.0f}/"
                   f"{r['cover']:.0f}/{r['dmm_cover']:.0f}" for p, by in rows.items()
                   for T, r in by.items()) + " ; " +
        " ; ".join(f"{p} base/cover/dmm_c {r['base']:.0f}/{r['cover']:.0f}/{r['dmm_cover']:.0f}"
                   for p, r in earlier.items() if "error" not in r))
    with open(os.path.join(args.out, "run.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
