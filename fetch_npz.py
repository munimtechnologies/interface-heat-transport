# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. (MV2). Made with AI agents working in MV2's name.
"""Download the tbc-tin-v1 orb-v3 npz files from the public repository, retrying until they are
published (the publish commit and this job may start in either order). Usage:
python fetch_npz.py <dest dir>"""
import os
import sys
import time
import urllib.request

BASE = ("https://raw.githubusercontent.com/munimtechnologies/interface-heat-transport/main/"
        "results/tin/orb-v3/")
FILES = [f"anh-orb-v3-{m}.npz" for m in ("TiN", "C", "Si", "Ge")]


def main():
    dest = sys.argv[1]
    os.makedirs(dest, exist_ok=True)
    for name in FILES:
        for attempt in range(90):
            try:
                urllib.request.urlretrieve(BASE + name, os.path.join(dest, name))
                with open(os.path.join(dest, name), "rb") as fh:
                    head = fh.read(2)
                if head != b"PK":
                    raise ValueError("not a zip/npz file")
                print(f"fetched {name} ({os.path.getsize(os.path.join(dest, name))} bytes) "
                      f"after {attempt} retries", flush=True)
                break
            except Exception as e:  # noqa: BLE001
                print(f"{name}: attempt {attempt}: {type(e).__name__}: {e}", flush=True)
                time.sleep(60)
        else:
            print(f"FAILED to fetch {name}", flush=True)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
