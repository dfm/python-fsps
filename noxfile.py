import os
from pathlib import Path

import nox

ALL_PYTHON_VS = ["3.8", "3.9", "3.10", "3.11"]


def _run_with_sps_home(session: nox.Session, *args, **kwargs):
    sps_home = os.environ.get(
        "SPS_HOME", Path(__file__).parent / "src" / "fsps" / "libfsps"
    )
    kwargs["env"] = dict(kwargs.get("env", {}), SPS_HOME=str(sps_home))
    return session.run(*args, **kwargs, external=True)


@nox.session(python=ALL_PYTHON_VS)
def tests(session):
    session.install(".[test]")
    _run_with_sps_home(session, "python", "tests/simple.py")
    _run_with_sps_home(
        session,
        "python",
        "-m",
        "pytest",
        "-n",
        "2",
        "--durations=0",
        "tests/tests.py",
    )


@nox.session(python=ALL_PYTHON_VS)
def options(session):
    session.install(
        ".[test]", env={"FFLAGS": "-DMIST=0 -DPADOVA=1 -DMILES=0 -DBASEL=1"}
    )
    _run_with_sps_home(session, "python", "tests/options.py")
