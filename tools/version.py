#!/usr/bin/env python3

from pathlib import Path

if Path(".git").exists():
    from setuptools_scm._cli import main

    main()

elif Path("src/fsps/fsps_version.py").exists():
    with open("src/fsps/fsps_version.py", "r") as f:
        exec(f.read())

    print(version) 

else:
    print("unknown")
