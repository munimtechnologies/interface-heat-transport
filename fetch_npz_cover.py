# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Download the orb-v3 npz files tbc_cover needs from the public repository into
<dest>/<repo path>. The four TiN-stage files are required; the earlier-stage files are optional
(a missing one is reported and its pair is skipped). Usage: python fetch_npz_cover.py <dest>"""
import os
import sys
import time
import urllib.request

BASE = "https://raw.githubusercontent.com/munimtechnologies/interface-heat-transport/main/"
REQUIRED = [f"results/tin/orb-v3/anh-orb-v3-{m}.npz" for m in ("TiN", "C", "Si", "Ge")]
OPTIONAL = [f"results/anh/v1/anh-orb-v3-{m}.npz" for m in ("Al", "Si", "GaN", "AlN")] + \
    [f"results/nonmetal/zno-v1/anh-orb-v3-{m}.npz" for m in ("ZnO", "GaN")]


def get(dest, path, tries):
    out = os.path.join(dest, path)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for attempt in range(tries):
        try:
            urllib.request.urlretrieve(BASE + path, out)
            with open(out, "rb") as fh:
                if fh.read(2) != b"PK":
                    raise ValueError("not a zip/npz file")
            print(f"fetched {path} ({os.path.getsize(out)} bytes)", flush=True)
            return True
        except Exception as e:  # noqa: BLE001
            print(f"{path}: attempt {attempt}: {type(e).__name__}: {e}", flush=True)
            time.sleep(20)
    if os.path.exists(out):
        os.remove(out)
    return False


def main():
    dest = sys.argv[1]
    for p in REQUIRED:
        if not get(dest, p, 10):
            print(f"FAILED to fetch required {p}", flush=True)
            return 1
    for p in OPTIONAL:
        if not get(dest, p, 5):
            print(f"optional {p} not fetched; its pair is skipped", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
