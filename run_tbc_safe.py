# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""run_tbc_safe: memory-safe driver for the stage-1 baselines of Grand Challenge 05 (job tbc-base-v2).

tbc-base-v1 was killed at its 3 GB cap after three of its ten materials (MACE-MPA-0 in float64 on
108-216-atom supercells). This driver computes each (potential, material) phonon spectrum in its own
spawned child process, so memory goes back to the system after every material and an out-of-memory
kill costs one material, not the whole job. If a MACE child is killed by a signal, that material is
retried once with MACE in float32 (recorded in info[name]['dtype']). Orb-v3 already runs in float32
and gets one attempt.

Physics, inputs, tables and grades are those of run_tbc.py, tbc_phonons.py and tbc_transport.py,
imported unchanged from the repository: this file changes no number. Last line of run.txt: 'TBC:'.
"""
import argparse
import json
import multiprocessing as mp
import os
import time
import traceback

import numpy as np

VERSION = "mv2 run-tbc-safe 0.1.0"

# Information only, not a pre-registered grade. Without LO-TO splitting the highest zone-centre mode
# of wurtzite GaN is the silent B1(high) mode, 85.5 meV = 20.67 THz (Ruf et al., Phys. Rev. Lett. 86,
# 906 (2001), inelastic x-ray scattering). PLAN.md's B2 reference for GaN (16.8 THz) is E1(TO).
REF_B1HIGH_THZ = {"GaN": 20.67}

INFO_KEYS = ("relax_converged", "cell_lengths", "n_atoms_per_m3", "max_thz", "min_thz",
             "imag_fraction", "seconds", "n_disp", "n_prim")


def attempts_for(model):
    """Precision attempts in order. MACE: float64 (as pre-registered), then float32 after a kill."""
    return ["float64", "float32"] if model == "mace-mpa-0" else ["default"]


def describe_exit(code):
    if code == 0:
        return "ok"
    if code == "timeout":
        return "timeout"
    if code is None:
        return "unknown"
    if code == -9:
        return "killed (SIGKILL, likely out of memory)"
    if code < 0:
        return f"signal {-code}"
    return f"exit {code}"


def should_retry(code):
    """Retry in lower precision only after a kill by a signal (memory), never after an exception."""
    return isinstance(code, int) and code < 0


def write_part(stem, arrays, info):
    np.savez(stem + ".npz", **arrays)
    tmp = stem + ".json.tmp"
    with open(tmp, "w") as fh:
        json.dump(info, fh)
    os.replace(tmp, stem + ".json")   # the .json appears only when the part is complete


def read_part(stem):
    with open(stem + ".json") as fh:
        info = json.load(fh)
    with np.load(stem + ".npz") as z:
        arrays = {k: np.array(z[k]) for k in z.files}
    return arrays, info


def finish_info(info, name, ref_max):
    info["ref_max_thz"] = ref_max
    info["max_ratio"] = (info["max_thz"] / ref_max) if ref_max else None
    if name in REF_B1HIGH_THZ:
        info["max_ratio_b1high_info"] = info["max_thz"] / REF_B1HIGH_THZ[name]
    return info


def load_calc(model, dtype):
    import tbc_phonons as tp
    if model == "mace-mpa-0" and dtype == "float32":
        from mace.calculators import mace_mp
        calc = mace_mp(model=tp.MACE_MPA0_URL, default_dtype="float32", device="cpu")
        print("MACE loaded in float32 (fallback after a kill)", flush=True)
        return calc, 0.01
    return tp.load_calculator(model)


def child(model, name, dtype, edges, normal, stem):
    """Runs in a fresh process: phonons of one material, its flux spectra and info -> stem.npz/.json."""
    import resource  # POSIX only; imported here so the module and its small tests load anywhere
    import tbc_phonons as tp
    import tbc_transport as tt
    calc, disp = load_calc(model, dtype)
    r = tp.material_phonons(name, calc, disp)
    arrays = {"iso": tt.flux_spectrum(r["freq"], r["vel"], r["weights"], r["vcell"], edges)}
    if normal is not None:
        arrays["or"] = tt.flux_spectrum(r["freq"], r["vel"], r["weights"], r["vcell"], edges,
                                        normal=normal)
    info = {k: r[k] for k in INFO_KEYS}
    info["relax_converged"] = bool(info["relax_converged"])
    info["n_disp"] = int(info["n_disp"])
    info["n_prim"] = int(info["n_prim"])
    info["dtype"] = dtype
    info["peak_rss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0, 3)
    write_part(stem, arrays, info)


def run_one(ctx, target, args, timeout_s):
    p = ctx.Process(target=target, args=args)
    p.start()
    p.join(max(1.0, timeout_s))
    if p.is_alive():
        p.terminate()
        p.join(10)
        if p.is_alive():
            p.kill()
            p.join()
        return "timeout"
    return p.exitcode


def run_model_safe(model, edges, out, deadline, log, ctx, mats, normals, refs, target=child):
    spectra, info = {}, {}
    work = os.path.join(out, "parts")
    os.makedirs(work, exist_ok=True)
    for name in mats:
        tried, done = [], False
        for dtype in attempts_for(model):
            left = deadline - time.time()
            if left < 30:
                tried.append(f"{dtype}:no time left")
                break
            stem = os.path.join(work, f"{model}-{name}-{dtype}")
            code = run_one(ctx, target, (model, name, dtype, edges, normals.get(name), stem), left)
            tried.append(f"{dtype}:{describe_exit(code)}")
            if code == 0 and os.path.exists(stem + ".json"):
                arrays, inf = read_part(stem)
                spectra[name] = arrays
                info[name] = finish_info(inf, name, refs.get(name))
                info[name]["attempts"] = list(tried)
                done = True
                break
            if not should_retry(code):
                break
        if not done:
            info[name] = {"failed": "; ".join(tried)}
            log(f"{model} {name}: FAILED {'; '.join(tried)}")
            continue
        d = info[name]
        log(f"{model} {name}: a,b,c={['%.4f' % x for x in d['cell_lengths']]} max={d['max_thz']:.2f} "
            f"THz (ref {d['ref_max_thz']}) min={d['min_thz']:.3f} imag={d['imag_fraction']:.4f} "
            f"conv={d['relax_converged']} {d['seconds']} s | {d['dtype']} peak {d['peak_rss_gb']} GB")
    np.savez_compressed(os.path.join(out, f"spectra-{model}.npz"), edges=edges,
                        **{f"{k}_{v}": a for k, dd in spectra.items() for v, a in dd.items()})
    return spectra, info


def safe_line(info_by_model):
    parts = []
    for m, info in info_by_model.items():
        ok = [n for n, d in info.items() if "failed" not in d]
        f32 = [n for n, d in info.items() if d.get("dtype") == "float32"]
        bad = [n for n, d in info.items() if "failed" in d]
        peak = max([d.get("peak_rss_gb", 0.0) for d in info.values() if "failed" not in d] or [0.0])
        g = info.get("GaN", {})
        b1 = g.get("max_ratio_b1high_info")
        parts.append(f"{m} done {len(ok)}/{len(info)} f32 [{','.join(f32)}] failed [{','.join(bad)}] "
                     f"peak {peak:.2f} GB GaN/B1high {('%.3f' % b1) if b1 else 'n/a'}")
    return "SAFE: " + " | ".join(parts)


def main():
    import run_tbc as rt
    import tbc_phonons as tp
    import tbc_transport as tt
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--models", default="mace-mpa-0,orb-v3")
    ap.add_argument("--budget-min", type=float, default=150.0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    deadline = t0 + args.budget_min * 60.0
    logf = open(os.path.join(args.out, "progress.txt"), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"{VERSION} | {rt.VERSION} | {tt.VERSION} | {tp.VERSION}")
    edges = tt.make_edges(48.0, 0.02)
    ctx = mp.get_context("spawn")
    res_by_model, info_by_model = {}, {}
    for m in [x.strip() for x in args.models.split(",") if x.strip()]:
        try:
            spectra, info = run_model_safe(m, edges, args.out, deadline, log, ctx, rt.MATS,
                                           rt.NORMALS, tp.REF_MAX_THZ)
        except Exception as e:
            log(f"MODEL {m} FAILED: {type(e).__name__}: {e}")
            log(traceback.format_exc())
            continue
        info_by_model[m] = info
        res_by_model[m] = rt.tables(spectra, info, edges)
    g = rt.grade(res_by_model, info_by_model)
    summary = {"version": rt.VERSION, "driver": VERSION, "transport": tt.VERSION, "grades": g,
               "info": info_by_model, "results": res_by_model,
               "eph": {"nominal": rt.EPH, "high": rt.EPH_HIGH},
               "b1high_reference_info": {"GaN_THz": REF_B1HIGH_THZ["GaN"],
                                         "source": "Ruf et al., PRL 86, 906 (2001), 85.5 meV"},
               "seconds": round(time.time() - t0, 1)}
    with open(os.path.join(args.out, "tbc-summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    for m, r in res_by_model.items():
        for kind in ("oriented", "iso"):
            for p, d in r[kind].items():
                row = " ".join(f"{T}:{d[T]['dmm']/1e6:.0f}/{d[T]['dmm_simons']/1e6:.0f}" for T in
                               ("100", "300", "400", "700"))
                log(f"G {m} {kind:8s} {p:8s} conv/simons MW/m2K {row} rad {d['300']['rad1']/1e6:.0f},"
                    f"{d['300']['rad2']/1e6:.0f}")
        for p, d in r["debye"].items():
            log(f"G {m} debye    {p:8s} 300K amm {d['300']['amm']/1e6:.0f} dmm {d['300']['dmm']/1e6:.0f}")
    for k, v in g.items():
        log(f"{k}: {v}")
    log(safe_line(info_by_model))
    log("TBC: " + " | ".join(f"{k} {v.split('(')[0]}" for k, v in g.items()))


if __name__ == "__main__":
    main()
