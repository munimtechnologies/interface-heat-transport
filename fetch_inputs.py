# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Fetch the unchanged stage-2a modules and tbc-anh-v1's saved orb-v3 lifetime files from the public
repository into /work. The lifetime files are committed right after the run that filed this job,
so missing files are retried for up to --wait-min minutes. Standard library only."""
import argparse
import os
import sys
import time
import urllib.request

BASE = "https://raw.githubusercontent.com/munimtechnologies/interface-heat-transport/main/"
MODULES = ("tbc_bte.py", "tbc_anh.py", "run_tbc_anh.py", "tbc_transport.py", "tbc_phonons.py",
           "test_tbc_bte.py", "test_tbc_transport.py")
MATS = ("Al", "Si", "GaN", "AlN", "C")
V1_DIR = "results/anh/v1/"


def get(url, dest):
    with urllib.request.urlopen(url, timeout=60) as resp:
        data = resp.read()
    with open(dest, "wb") as fh:
        fh.write(data)
    return len(data)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/work")
    ap.add_argument("--wait-min", type=float, default=90.0)
    ap.add_argument("--modules", action="store_true", help="also fetch the stage-2a modules")
    args = ap.parse_args()
    if args.modules:
        for p in MODULES:
            print(f"{p}: {get(BASE + p, os.path.join(args.work, p))} bytes", flush=True)
    dest = os.path.join(args.work, "v1parts")
    os.makedirs(dest, exist_ok=True)
    todo = [f"anh-orb-v3-{m}.npz" for m in MATS]
    deadline = time.time() + args.wait_min * 60.0
    while todo:
        left = []
        for f in todo:
            try:
                n = get(BASE + V1_DIR + f, os.path.join(dest, f))
                print(f"{f}: {n} bytes", flush=True)
            except Exception as e:  # not committed yet, or a network error: retry
                left.append(f)
                print(f"{f}: not yet ({type(e).__name__})", flush=True)
        todo = left
        if todo:
            if time.time() > deadline:
                print("MISSING after waiting: " + ", ".join(todo), flush=True)
                return 1
            time.sleep(120)
    return 0


if __name__ == "__main__":
    sys.exit(main())
