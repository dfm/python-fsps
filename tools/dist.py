#!/usr/bin/env python3

import glob
import os
from pathlib import Path
from shutil import copyfile, rmtree

dist = Path(os.environ["MESON_DIST_ROOT"])
source = Path(os.environ["MESON_SOURCE_ROOT"])

print("Removing FSPS data files from dist...")
for path in glob.glob(str(dist / "src" / "fsps" / "libfsps" / "*")):
    path = Path(path)
    if path.name in ("src", "LICENSE", "README.md"):
        continue
    print(f"- {path.name}")
    if path.is_file():
        path.unlink()
    else:
        rmtree(path)

print("Adding version file to dist...")
copyfile(
    source / "src" / "fsps" / "fsps_version.py",
    dist / "src" / "fsps" / "fsps_version.py",
)
