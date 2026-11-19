# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""tbc_tin: stage 3b of the interface-heat-transport study (PLAN-amendment-0.5.2.md).

Khan, Shi, Feser and Wilson (Appl. Phys. Lett. 125, 2024; OSTI 2578597) measured the thermal
boundary conductance of epitaxial TiN on diamond, SiC, Si and Ge from about 120 to 1000 K. It
"increases monotonically without stagnating at higher temperatures". These are the highest-
temperature interface data we know of, and the hardest test of a no-fit model.

This module holds
  * the blind reading of that paper (and of Koh et al., PRB 102, 205304, 2020) as constants,
    checked against data/read-hight.json by the tests;
  * exploratory data-only numbers (numbers()): how far a purely elastic conductance can rise,
    conventional (model-free heat-capacity bound) and dressed (Simons / Shi gray formula with the
    paper's own curves);
  * the registration of TiN (rocksalt) and Ge (diamond) and of (001) interface geometries in the
    stage 2a/2b tables (register_tin), so run_tbc_tin.py can reuse that code unchanged;
  * grade_t: the pre-registered grades T0-T6 of job tbc-tin-v1, written before it runs.

Pure Python apart from register_tin (which imports the stage 2 modules only when called).
"""
import math

VERSION = "mv2 tbc-tin 0.1.0"

H_OVER_K_K_PER_THZ = 47.99243  # h / k_B in K per THz

# ---- blind reading (MV2's hands, 2026-11-19; data/read-hight.json) --------------------------
# Khan et al. 2024, measured G (MW m^-2 K^-1) at 300 ... K, linear interpolation between points.
KHAN = {
    "TiN/C": {300: 449.0, 400: 572.0, 500: 662.0, 600: 740.0, 700: 814.0, 800: 870.0, 900: 910.0},
    "TiN/SiC": {300: 748.0, 400: 842.0, 500: 900.0, 600: 920.0, 700: 968.0},
    "TiN/Si": {300: 507.0, 400: 633.0, 500: 736.0, 600: 870.0, 700: 916.0},
    "TiN/Ge": {300: 190.0, 400: 228.0, 500: 244.0, 600: 263.0, 700: 281.0},
}
KHAN_LAST = {"TiN/C": (996.0, 948.0), "TiN/SiC": (771.0, 999.0), "TiN/Si": (771.0, 942.0),
             "TiN/Ge": (771.0, 299.0)}
# model curves drawn in the paper, values at (300 K, 1000 K)
KHAN_CURVES = {
    "TiN/C": {"elastic_limit": (364.0, 505.0), "DMM": (259.0, 352.0),
              "R_TiN": (1290.0, 1510.0), "R_other": (3200.0, 8200.0)},
    "TiN/Si": {"elastic_limit": (411.0, 438.0), "DMM": (291.0, 311.0),
               "R_TiN": (1280.0, 1500.0), "R_other": (1065.0, 1230.0)},
}
# reading uncertainty: +-5-15 on clear markers, +-20-30 on the partly hidden Si points 480-771 K
READ_ERR = {"TiN/C": 15.0, "TiN/SiC": 15.0, "TiN/Si": 30.0, "TiN/Ge": 15.0}
# Koh et al. 2020 (vector figure, +-1): measured at 300 / 400 / ~500 K and the modal
# non-equilibrium Landauer curve with DMM transmission (Shi's method) at 300 / 500 K
KOH = {"Al/sapphire": (315.0, 338.0, 325.0), "Co/sapphire": (370.0, 405.0, 460.0),
       "Ru/sapphire": (351.0, 340.0, 361.0)}
KOH_NEL = {"Al/sapphire": (334.0, 351.0), "Co/sapphire": (372.0, 390.0),
           "Ru/sapphire": (296.0, 307.0)}

# highest phonon frequency (THz) used for the model-free bound; reference values
NU_MAX_REF = {"TiN": 18.0, "Si": 15.6, "Ge": 9.1, "C": 39.9}

# the temperature each pair is graded at (highest 100 K step inside the measured range)
T_HIGH = {"TiN/C": 900, "TiN/Si": 700, "TiN/Ge": 700}
PAIRS = ("TiN/C", "TiN/Si", "TiN/Ge")


# ---- small physics helpers ------------------------------------------------------------------
def cap(nu_thz, T):
    """Einstein heat capacity of one mode, in units of k_B."""
    x = H_OVER_K_K_PER_THZ * nu_thz / T
    if x < 1e-8:
        return 1.0
    ex = math.exp(x)
    return x * x * ex / (ex - 1.0) ** 2


def cap_ratio(nu_thz, T1, T2):
    return cap(nu_thz, T2) / cap(nu_thz, T1)


def dressed_gray(g0, r1, r2):
    """Simons / Zeng-Chen / Shi gray dressed conductance, 1/G_d = 1/G_0 - 1/2R_1 - 1/2R_2.
    Holds for any gray transmission obeying detailed balance (G_0 = t12 R1 = t21 R2)."""
    inv = 1.0 / g0 - 0.5 / r1 - 0.5 / r2
    return float("inf") if inv <= 0 else 1.0 / inv


def floor_nonelastic(g1, g2, kappa):
    """Least non-elastic conductance at T2 if the elastic part grows at most kappa-fold from T1
    (and the non-elastic part at T1 is >= 0): G2 - kappa G1."""
    return g2 - kappa * g1


def channel_floor_linear(g1, g2, rho, kappa):
    """Least T-proportional channel at T1 (I = c T, rho = T2/T1) if E2 <= kappa E1."""
    return (g2 - kappa * g1) / (rho - kappa)


def interp(table, T):
    ks = sorted(table)
    if T in table:
        return table[T]
    for a, b in zip(ks, ks[1:]):
        if a < T < b:
            return table[a] + (table[b] - table[a]) * (T - a) / (b - a)
    raise ValueError(T)


# ---- exploratory data-only numbers (amendment 0.5.2, section 2) ------------------------------
def numbers():
    x = {}
    # Khan's own curves, dressed (gray), at 300 and 1000 K
    for pair, c in KHAN_CURVES.items():
        for name in ("elastic_limit", "DMM"):
            for i, T in enumerate((300, 1000)):
                x[f"{pair} dressed {name} {T}"] = dressed_gray(c[name][i], c["R_TiN"][i],
                                                               c["R_other"][i])
    gC1, (TC2, gC2) = KHAN["TiN/C"][300], KHAN_LAST["TiN/C"]
    gS1, (TS2, gS2) = KHAN["TiN/Si"][300], KHAN_LAST["TiN/Si"]
    kC = x["TiN/C dressed elastic_limit 1000"] / x["TiN/C dressed elastic_limit 300"]
    kS = x["TiN/Si dressed elastic_limit 1000"] / x["TiN/Si dressed elastic_limit 300"]
    x["TiN/C kappa dressed EL 300-1000"] = kC
    x["TiN/Si kappa dressed EL 300-1000"] = kS
    x["TiN/C excess over dressed EL at 996"] = gC2 - x["TiN/C dressed elastic_limit 1000"]
    x["TiN/C floor nonelastic 996 (dressed EL trend)"] = floor_nonelastic(gC1, gC2, kC)
    x["TiN/Si floor nonelastic 771 (dressed EL trend, 300-1000 ratio)"] = floor_nonelastic(gS1, gS2, kS)
    x["TiN/C linear channel floor I300 (dressed EL trend)"] = channel_floor_linear(
        gC1, gC2, TC2 / 300.0, kC)
    # model-free conventional elastic bound: growth <= heat-capacity ratio of the highest mode
    # that can cross elastically (the softer side's top frequency)
    for pair, side, T2, g2 in (("TiN/C", "TiN", TC2, gC2), ("TiN/Si", "Si", TS2, gS2)):
        k = cap_ratio(NU_MAX_REF[side], 300.0, T2)
        x[f"{pair} kappa conv max ({side} {NU_MAX_REF[side]} THz) 300-{int(T2)}"] = k
        x[f"{pair} floor nonelastic {int(T2)} (model-free conventional)"] = floor_nonelastic(
            KHAN[pair][300], g2, k)
    kG = cap_ratio(NU_MAX_REF["Ge"], 300.0, 771.0)
    x["TiN/Ge kappa conv max (Ge 9.1 THz) 300-771"] = kG
    x["TiN/Ge floor nonelastic 771 (model-free conventional)"] = floor_nonelastic(
        KHAN["TiN/Ge"][300], KHAN_LAST["TiN/Ge"][1], kG)
    # measured rises used by T1/T3
    for p in PAIRS:
        x[f"{p} measured ratio {T_HIGH[p]}/300"] = KHAN[p][T_HIGH[p]] / KHAN[p][300]
    return x


# ---- registration for run_tbc_tin.py (stage 2a/2b code reused, nothing copied) ---------------
TIN = dict(kind="rocksalt", el="TiN", a=4.2417, sc=(3, 3, 3), mesh=28)
GE = dict(kind="diamond", el="Ge", a=5.6579, sc=(3, 3, 3), mesh=28)
MATS = ["TiN", "C", "Si", "Ge"]
PAIR_TUPLES = [("TiN", "C"), ("TiN", "Si"), ("TiN", "Ge")]
TEMPS = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0]
G_EP_TIN = 1.0e18   # W m^-3 K^-1, order of magnitude quoted by Khan et al.; sensitivity only

# terminating (001) layers, same convention as tbc_inel.GEOMETRY:
#   rocksalt(001): 6 neighbours, 1 across; sigma = 4/a^2, n = 8/a^3 -> sigma/n = a/2
#   diamond(001):  4 bonds, 2 across;      sigma = 2/a^2, n = 8/a^3 -> sigma/n = a/4
GEOM_001 = {"rocksalt001": (1.0 / 6.0, 0.5), "diamond001": (2.0 / 4.0, 0.25)}


def t_eff_001(geometry, a):
    f, s = GEOM_001[geometry]
    return f * s * a


def register_tin():
    import tbc_anh as ta
    import tbc_inel as ti
    import tbc_phonons as tp
    import run_tbc_anh as ra
    import run_tbc_inel as ri
    tp.MATERIALS["TiN"] = dict(TIN)
    tp.MATERIALS["Ge"] = dict(GE)
    tp.REF_MAX_THZ["TiN"] = NU_MAX_REF["TiN"]
    tp.REF_MAX_THZ["Ge"] = NU_MAX_REF["Ge"]
    ta.ANH["TiN"] = dict(sc3=(2, 2, 2), cut=4.0, mesh=(16, 16, 16), lattice="cubic")
    ta.ANH["Ge"] = dict(sc3=(2, 2, 2), cut=None, mesh=(16, 16, 16), lattice="cubic")
    ta.KAPPA_RANGE_300["TiN"] = (1.0, 300.0)
    ta.KAPPA_RANGE_300["Ge"] = (25.0, 120.0)
    for name, (f, s) in GEOM_001.items():
        ti.GEOMETRY[name] = dict(f=f, sigma_over_n=(lambda a, c, s=s: s * a))
    for m in MATS:
        ra.NORMALS[m] = [0, 0, 1]
        ra.KDIR[m] = "iso"
        ri.GEOM[m] = "rocksalt001" if m == "TiN" else "diamond001"
        ri.NCONV[m] = 8
    ra.MATS = list(MATS)
    ra.PAIRS = list(PAIR_TUPLES)
    ra.TEMPS = list(TEMPS)


# ---- pre-registered grades of tbc-tin-v1 (amendment 0.5.2, section 3) ------------------------
T0_MAX_THZ = {"TiN": (14.0, 23.0), "C": (35.0, 45.0), "Si": (13.5, 17.5), "Ge": (7.5, 10.5)}
T0_IMAG_MAX = 0.02
T2_TOL = 0.25
T3_TOL = 0.15
T6_RANGE = (0.25, 0.55)


def _rows(prim, pair, T):
    return prim.get("pairs", {}).get(pair, {}).get(str(int(T)))


def grade_t(prim):
    """prim = tbt-summary.json['results'][primary model]. Returns {T0..T6: 'PASS(...)'...}."""
    g = {}
    mats = prim.get("materials", {})
    bad = []
    for m, (lo, hi) in T0_MAX_THZ.items():
        v = mats.get(m)
        if v is None:
            bad.append(f"{m} missing")
            continue
        if not (lo <= v["max_thz"] <= hi):
            bad.append(f"{m} max {v['max_thz']:.2f}")
        if v["imag_fraction"] > T0_IMAG_MAX:
            bad.append(f"{m} imag {v['imag_fraction']:.3f}")
    for p in PAIRS:
        if _rows(prim, p, 300) is None or _rows(prim, p, T_HIGH[p]) is None:
            bad.append(f"{p} rows missing")
    g["T0"] = "PASS" if not bad else "FAIL(" + "; ".join(bad) + ")"

    have = [p for p in PAIRS if _rows(prim, p, 300) and _rows(prim, p, T_HIGH[p])]
    if len(have) < len(PAIRS):
        for k in ("T1", "T2", "T3", "T4", "T5", "T6"):
            g[k] = "UNDETERMINED(pairs missing: " + ",".join(sorted(set(PAIRS) - set(have))) + ")"
        return g
    lo = {p: _rows(prim, p, 300) for p in PAIRS}
    hi = {p: _rows(prim, p, T_HIGH[p]) for p in PAIRS}
    meas_r = {p: KHAN[p][T_HIGH[p]] / KHAN[p][300] for p in PAIRS}

    # T1: the model's elastic (dressed BTE) conductance rises less than the data, every pair
    el_r = {p: hi[p]["G_bte"] / lo[p]["G_bte"] for p in PAIRS}
    ok = all(el_r[p] < meas_r[p] for p in PAIRS)
    g["T1"] = ("PASS" if ok else "FAIL") + "(" + ", ".join(
        f"{p} {el_r[p]:.3f} vs {meas_r[p]:.3f}" for p in PAIRS) + ")"

    # T2: no-fit total within 25 % of the data at 300 K and at T_high, for >= 2 of 3 pairs
    n2, parts = 0, []
    for p in PAIRS:
        d1 = lo[p]["G_total"] / KHAN[p][300] - 1.0
        d2 = hi[p]["G_total"] / KHAN[p][T_HIGH[p]] - 1.0
        n2 += int(abs(d1) <= T2_TOL and abs(d2) <= T2_TOL)
        parts.append(f"{p} {100 * d1:+.0f}%/{100 * d2:+.0f}%")
    g["T2"] = ("PASS" if n2 >= 2 else "FAIL") + f"({n2}/3: " + ", ".join(parts) + ")"

    # T3: the no-fit total rises like the data (ratio within 15 %), for >= 2 of 3 pairs
    n3, parts = 0, []
    for p in PAIRS:
        r = hi[p]["G_total"] / lo[p]["G_total"]
        n3 += int(abs(r / meas_r[p] - 1.0) <= T3_TOL)
        parts.append(f"{p} {r:.3f} vs {meas_r[p]:.3f}")
    g["T3"] = ("PASS" if n3 >= 2 else "FAIL") + f"({n3}/3: " + ", ".join(parts) + ")"

    # T4: closer to the data than the model's own DMM at T_high, every pair
    ok, parts = True, []
    for p in PAIRS:
        m = KHAN[p][T_HIGH[p]]
        e_tot, e_dmm = abs(hi[p]["G_total"] - m), abs(hi[p]["G_dmm"] - m)
        ok &= e_tot < e_dmm
        parts.append(f"{p} {e_tot:.0f} vs {e_dmm:.0f}")
    g["T4"] = ("PASS" if ok else "FAIL") + "(" + ", ".join(parts) + ")"

    # T5: ordering at 300 K, Si > C > Ge as measured (507 > 449 > 190)
    a, b, c = (lo[p]["G_total"] for p in ("TiN/Si", "TiN/C", "TiN/Ge"))
    g["T5"] = ("PASS" if a > b > c else "FAIL") + f"(Si {a:.0f}, C {b:.0f}, Ge {c:.0f})"

    # T6: Ge / Si at 300 K in [0.25, 0.55] (measured 0.375)
    q = lo["TiN/Ge"]["G_total"] / lo["TiN/Si"]["G_total"]
    g["T6"] = ("PASS" if T6_RANGE[0] <= q <= T6_RANGE[1] else "FAIL") + f"({q:.3f})"
    return g
