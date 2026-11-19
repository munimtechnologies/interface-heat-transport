# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""run_anh_safe: memory-safe driver for stage 2a (job tbc-anh-v2, PLAN-amendment-0.2.md and 0.5.5).

tbc-anh-v1 finished the orb-v3 lifetimes of all five materials, then was killed at its 3 GB cap
in the MACE-MPA-0 cross-check (third material, GaN). This driver computes each (potential,
material) in its own spawned child process, so memory goes back to the system after every
material and a kill costs one material. A MACE child killed by a signal is retried once in float32
(recorded). Each child saves anh-<model>-<material>.npz with the same arrays as v1 plus the cell
volume, and a .json that appears only when the part is complete.

Physics and grades are those of tbc_anh.py and run_tbc_anh.py, imported unchanged; the tables and
grades are made by anh_pairs.compute from the saved files (orb-v3 from v1, MACE from this job).
Last line of run.txt starts with 'TBA:'.
"""
import argparse
import json
import multiprocessing as mp
import os
import time

import numpy as np

VERSION = "mv2 run-anh-safe 0.1.0"


def attempts_for(model):
    """Precision attempts in order. MACE: float64 (as in v1), then float32 only after a kill."""
    return ["float64", "float32"] if model == "mace-mpa-0" else ["default"]


def should_retry(code):
    """Retry in lower precision only after a kill by a signal (memory), never after an exception."""
    return isinstance(code, int) and code < 0


def describe_exit(code):
    if code == 0:
        return "ok"
    if code == "timeout":
        return "timeout"
    if code is None:
        return "unknown"
    if code == -9:
        return "killed (SIGKILL, likely out of memory)"
    if isinstance(code, int) and code < 0:
        return f"signal {-code}"
    return f"exit {code}"


def child(model, name, dtype, stem):
    """Runs in a fresh process: lifetimes and flux spectra of one material -> stem.npz + stem.json."""
    import resource  # POSIX only; imported here so the module and its tests load anywhere
    import run_tbc_anh as rta
    import tbc_anh as ta
    import tbc_phonons as tp
    import tbc_transport as tt
    if model == "mace-mpa-0" and dtype == "float32":
        from mace.calculators import mace_mp
        calc = mace_mp(model=tp.MACE_MPA0_URL, default_dtype="float32", device="cpu")
        print("MACE loaded in float32 (retry after a kill)", flush=True)
    else:
        calc, _ = tp.load_calculator(model)
    edges_f = tt.make_edges(48.0, 0.02)
    edges_c = tt.make_edges(48.0, 0.2)
    r = ta.anharmonic(name, calc, rta.TEMPS, log=lambda s: print(s, flush=True))
    fm = ta.flux_mesh(r)
    M_or = tt.flux_spectrum(fm["freq"], fm["vel"], fm["weights"], fm["vcell"], edges_f,
                            normal=rta.NORMALS[name])
    M_iso = tt.flux_spectrum(fm["freq"], fm["vel"], fm["weights"], fm["vcell"], edges_f)
    np.savez_compressed(stem + ".npz", edges_c=edges_c, freq=r["freq"], weights=r["weights"],
                        gamma3=r["gamma3"], gamma_iso=r["gamma_iso"], v2_iso=r["v2_iso"],
                        v2_z=r["v2_z"], kappa_ph3=r["kappa_ph3"], temps=np.array(rta.TEMPS),
                        M_or=M_or, M_iso=M_iso, edges_f=edges_f, vcell=np.array(r["vcell"]),
                        cell_lengths=np.array(r["cell_lengths"]))
    info = {"model": model, "material": name, "dtype": dtype,
            "relax_converged": bool(r["relax_converged"]), "cell_lengths": r["cell_lengths"],
            "n_fc3": int(r["n_fc3"]), "seconds": r["seconds"],
            "seconds_forces": r["seconds_forces"],
            "min_thz_dense": float(fm["freq"].min()), "max_thz_dense": float(fm["freq"].max()),
            "imag_fraction_dense": float(np.mean(fm["freq"] < -0.05)),
            "peak_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0, 3)}
    tmp = stem + ".json.tmp"
    with open(tmp, "w") as fh:
        json.dump(info, fh)
    os.replace(tmp, stem + ".json")


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


def run_model_safe(model, mats, parts, deadline, log, ctx, target=child):
    """Each material in its own child; returns {material: info or {'failed': ...}}."""
    os.makedirs(parts, exist_ok=True)
    info = {}
    for name in mats:
        tried, done = [], False
        stem = os.path.join(parts, f"anh-{model}-{name}")
        for dtype in attempts_for(model):
            left = deadline - time.time()
            if left < 30:
                tried.append(f"{dtype}:no time left")
                break
            code = run_one(ctx, target, (model, name, dtype, stem), left)
            tried.append(f"{dtype}:{describe_exit(code)}")
            if code == 0 and os.path.exists(stem + ".json"):
                with open(stem + ".json") as fh:
                    info[name] = json.load(fh)
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
        log(f"{model} {name}: done ({d['dtype']}) a,b,c={['%.4f' % x for x in d['cell_lengths']]} "
            f"fc3 {d['n_fc3']} {d['seconds']} s peak {d['peak_rss_gb']} GB imag {d['imag_fraction_dense']:.4f}")
    return info


def main():
    import anh_pairs as apr
    import run_tbc_anh as rta
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--models", default="mace-mpa-0", help="models to compute here")
    ap.add_argument("--reuse", action="append", default=[],
                    help="MODEL=DIR of saved files to grade first (orb-v3 from tbc-anh-v1)")
    ap.add_argument("--budget-min", type=float, default=170.0)
    ap.add_argument("--nodes", type=int, default=400)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    deadline = t0 + args.budget_min * 60.0
    logf = open(os.path.join(args.out, "progress.txt"), "w")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()

    log(f"{VERSION} | {apr.VERSION} | {rta.VERSION}")
    ctx = mp.get_context("spawn")
    parts = os.path.join(args.out, "parts")
    info = {}
    for m in [x.strip() for x in args.models.split(",") if x.strip()]:
        info[m] = run_model_safe(m, rta.MATS, parts, deadline, log, ctx)
    with open(os.path.join(args.out, "safe-info.json"), "w") as fh:
        json.dump(info, fh, indent=1)
    sources = apr.parse_sources(args.reuse) + [(m, parts) for m in info]
    ok = {m: [n for n, d in v.items() if "failed" not in d] for m, v in info.items()}
    log("SAFE: " + " | ".join(f"{m} done {len(v)}/{len(rta.MATS)} [{','.join(v)}]" for m, v in ok.items()))
    # tables and grades for all sources, written by anh_pairs (same output format as run_tbc_anh)
    argv = []
    for m, folder in sources:
        argv += ["--source", f"{m}={folder}"]
    argv += ["--out", args.out, "--nodes", str(args.nodes)]
    apr.main(argv)


if __name__ == "__main__":
    main()
