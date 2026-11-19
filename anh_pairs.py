# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""anh_pairs: stage-2a interface tables and grades (PLAN-amendment-0.2.md) from saved per-material
lifetime files anh-<model>-<material>.npz, without recomputing any force constant.

Job tbc-anh-v1 computed and saved the orb-v3 lifetimes of all five materials (Al, Si, GaN, AlN,
diamond) and was then killed at its 3 GB memory cap during the MACE-MPA-0 cross-check, before it
wrote its grades. This module rebuilds exactly what run_tbc_anh.py would have done next: the
per-bin conductivities (tbc_anh.kappa_bins), the interface tables (run_tbc_anh.pair_tables) and the
grades A0, A1, S2, S3, S5, S6 (run_tbc_anh.grade), all imported unchanged. It changes no number.

One input is not in the v1 files: the primitive-cell volume. For files that lack it we use the
relaxed cell lengths that tbc-anh-v1 printed (progress.txt, 4 decimals). A0 (per-bin sum against
phono3py's own conductivity, which used the exact volume) checks that choice.

Usage: python anh_pairs.py --source orb-v3=DIR [--source mace-mpa-0=DIR] --out OUT
Last line of run.txt starts with 'TBA:'.
"""
import argparse
import hashlib
import json
import os
import time

import numpy as np

import run_tbc_anh as rta
import tbc_anh as ta
import tbc_transport as tt

VERSION = "mv2 anh-pairs 0.1.0"

# Relaxed cell lengths (Angstrom) printed by tbc-anh-v1 (out/progress.txt) for the files that do not
# store the cell volume: conventional cubic cell for Al, Si, C; hexagonal cell for GaN, AlN.
V1_CELLS = {
    "orb-v3": {"Al": (4.0245, 4.0245, 4.0245), "Si": (5.4608, 5.4608, 5.4608),
               "GaN": (3.2180, 3.2180, 5.2454), "AlN": (3.1320, 3.1320, 5.0086),
               "C": (3.5721, 3.5721, 3.5721)},
    "mace-mpa-0": {"Al": (4.0387, 4.0387, 4.0387), "Si": (5.4674, 5.4674, 5.4674),
                   "GaN": (3.2047, 3.2047, 5.2503)},
}
HEX = ("GaN", "AlN")


def prim_volume(name, lengths):
    """Primitive-cell volume (m^3): a^3/4 for the fcc and diamond conventional cubic cell,
    (sqrt 3 / 2) a^2 c for the wurtzite hexagonal cell (which is primitive)."""
    a, b, c = (float(x) for x in lengths)
    if name in HEX:
        return 0.5 * np.sqrt(3.0) * a * b * c * 1e-30
    return a * b * c / 4.0 * 1e-30


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_material(path, model, name):
    """Arrays of one saved material, with r['vcell'] (m^3) and the source of the volume."""
    with np.load(path) as z:
        r = {k: np.array(z[k]) for k in z.files}
    if "vcell" in r:
        r["vcell"] = float(r["vcell"])
        src = "file"
    else:
        r["vcell"] = prim_volume(name, V1_CELLS[model][name])
        src = "printed cell lengths"
    return r, src


def material_entry(r, name, edges_c):
    """The per-material structure run_tbc_anh.per_material builds: kappa sums and per-bin arrays
    for every temperature and the needed directions, plus the flux spectra."""
    d = {"kappa": {}, "bins": {}, "M_or": r["M_or"], "M_iso": r["M_iso"]}
    for it, T in enumerate(float(t) for t in r["temps"]):
        key = str(int(T))
        kd = {}
        for dirn in sorted({"iso", rta.KDIR[name]}):
            k_iso = ta.kappa_bins(r, it, T, edges_c, dirn, with_isotope=True)
            k_3 = ta.kappa_bins(r, it, T, edges_c, dirn, with_isotope=False)
            k_cl = ta.kappa_bins(r, it, T, edges_c, dirn, with_isotope=False, classical=True)
            d["bins"][(key, dirn)] = (k_iso, k_3, k_cl)
            kd[dirn] = {"sum_bins_iso": float(k_iso.sum()), "sum_bins_3ph": float(k_3.sum()),
                        "phono3py_iso": ta.kappa_ph3_dir(r, it, dirn)}
        d["kappa"][key] = kd
    return d


def compute(sources, nodes=400, log=print):
    """sources: list of (model, folder). Returns per_model, tables, models, grades, file info."""
    edges_f = tt.make_edges(48.0, 0.02)
    edges_c = tt.make_edges(48.0, 0.2)
    per_model, tables, models, files = {}, {}, [], {}
    for m, folder in sources:
        data = {}
        for name in rta.MATS:
            path = os.path.join(folder, f"anh-{m}-{name}.npz")
            if not os.path.exists(path):
                log(f"{m} {name}: no file {path}")
                continue
            r, src = load_material(path, m, name)
            if r["edges_c"].shape != edges_c.shape or not np.allclose(r["edges_c"], edges_c) \
                    or r["edges_f"].shape != edges_f.shape or not np.allclose(r["edges_f"], edges_f):
                log(f"{m} {name}: frequency bins differ from tbc_transport.make_edges, skipped")
                continue
            d = material_entry(r, name, edges_c)
            data[name] = d
            files[f"{m}/{name}"] = {"sha256": sha256(path), "vcell_m3": r["vcell"],
                                    "vcell_source": src}
            k300 = d["kappa"]["300"]
            log(f"{m} {name}: kappa300 " + ", ".join(
                f"{dn} bins {v['sum_bins_iso']:.1f} (3ph only {v['sum_bins_3ph']:.1f}) phono3py "
                f"{v['phono3py_iso']:.1f}" for dn, v in k300.items()) + f" W/mK; volume from {src}")
        if not data:
            continue
        models.append(m)
        per_model[m] = data
        tables[m] = rta.pair_tables(data, edges_f, edges_c, nodes)
    g = rta.grade(per_model, tables, models)
    return per_model, tables, models, g, files


def parse_sources(items):
    out = []
    for s in items:
        m, _, folder = s.partition("=")
        if not m or not folder:
            raise ValueError(f"--source needs MODEL=DIR, got {s!r}")
        out.append((m.strip(), folder.strip()))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", action="append", required=True, help="MODEL=DIR, in grading order")
    ap.add_argument("--out", required=True)
    ap.add_argument("--nodes", type=int, default=400)
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    logf = open(os.path.join(args.out, "pairs-progress.txt"), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    import tbc_bte as tb
    log(f"{VERSION} | {rta.VERSION} | {tb.VERSION} | {ta.VERSION} | {tt.VERSION}")
    per_model, tables, models, g, files = compute(parse_sources(args.source), args.nodes, log)
    gep = tt.g_ep(*rta.EPH_AL)
    summary = {"version": rta.VERSION, "driver": VERSION, "bte": tb.VERSION, "grades": g,
               "models": models, "files": files,
               "info": {m: {n: {"kappa": d["kappa"]} for n, d in per_model[m].items()}
                        for m in models},
               "tables": tables, "g_ep_Al": gep, "seconds": round(time.time() - t0, 1)}
    with open(os.path.join(args.out, "tba-summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    for m in models:
        for p, rows in tables[m].items():
            line = " ".join(
                f"{T}:{rows[T]['iso']['G_dmm']/1e6:.0f}/{rows[T]['iso']['G']/1e6:.0f}/"
                f"{rows[T]['3ph']['G']/1e6:.0f}/{rows[T]['classical']['G']/1e6:.0f}/"
                f"{rows[T]['iso']['G_simons']/1e6:.0f}"
                for T in ("100", "300", "400", "500", "700") if T in rows)
            log(f"G {m} {p:18s} T:dmm/bte/bte3ph/bteclassical/simons MW/m2K {line}")
            if p.startswith("oriented Al/") and "300" in rows:
                s = tt.series(rows["300"]["iso"]["G"], gep)
                log(f"  {m} {p} 300 K with Al e-ph in series: {s/1e6:.0f} MW/m2K "
                    f"(measured {rta.MEAS[p.split()[1]]/1e6:.0f})")
    for k, v in g.items():
        log(f"{k}: {v}")
    log("TBA: " + " | ".join(f"{k} {v.split('(')[0]}" for k, v in g.items()))


if __name__ == "__main__":
    main()
