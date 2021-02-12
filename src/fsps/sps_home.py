# -*- coding: utf-8 -*-

__all__ = ["check_sps_home"]

import os


def check_sps_home():
    if "SPS_HOME" not in os.environ:
        raise RuntimeError("You need to have the SPS_HOME environment variable")

    path = os.environ["SPS_HOME"]
    if not os.path.isdir(path):
        raise RuntimeError(f"SPS_HOME environment variable '{path}' is not a directory")

    if not os.path.exists(os.path.join(path, "data", "emlines_info.dat")):
        raise RuntimeError(
            f"The FSPS directory at {path} doesn't seem to have the right data "
            "files installed"
        )
