#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""feature_demo.py - This script creates a set of PDFs that illustrate the
effect on the SED of successively turning on various options or changing the
value of some variables.
"""

import os, time

import matplotlib.pyplot as pl
from matplotlib.backends.backend_pdf import PdfPages

import fsps


def makefig(sps, tage=13.7, oldspec=None, **plotkwargs):
    w, spec = sps.get_spectrum(tage=tage)
    fig, ax = pl.subplots()
    if oldspec is not None:
        ax.plot(w, oldspec / w * 1e19, color="gray", linewidth=2, alpha=0.5)
    ax.plot(w, spec / w * 1e19, "C2", linewidth=2)
    return fig, ax, spec


def prettify(fig, ax, label=None):
    ax.set_xlim(0.9e3, 1e6)
    ax.set_xscale("log")
    ax.set_ylim(0.01, 2)
    # ax.set_yscale('log')
    ax.set_xlabel(r"rest-frame $\lambda$ ($\AA$)", fontsize=20)
    ax.set_ylabel(r"$\lambda \, f_\lambda$", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=16)
    if label is not None:
        ax.text(0.63, 0.85, label, transform=ax.transAxes, fontsize=16)

    fig.tight_layout()
    return fig, ax


def get_githash(dir):
    """
    Get the current git hash for a repo
    """
    import subprocess
    cmd = ["git", "rev-parse", "HEAD"]
    hash = subprocess.check_output(cmd, cwd=dir).decode("utf-8").strip()
    return hash


if __name__ == "__main__":

    pl.rcParams["font.family"] = "serif"
    pl.rcParams["font.serif"] = ["cmr10"]
    pl.rcParams["mathtext.fontset"] = "cm"
    pl.rcParams["axes.formatter.use_mathtext"] = True
    pl.rcParams["axes.grid"] = True
    pl.rcParams["grid.alpha"] = 0.5
    pl.rcParams["xtick.direction"] = "in"
    pl.rcParams["ytick.direction"] = "in"


    pyfsps_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    fsps_hash = get_githash(os.environ["SPS_HOME"])
    pyfsps_hash = get_githash(os.path.dirname(__file__))
    pyfsps_lib_hash = get_githash(os.path.join(pyfsps_dir, "src", "fsps", "libfsps"))
    print(f"FSPS git hash: {fsps_hash}")
    print(f"pyFSPS git hash: {pyfsps_hash}")
    print(f"pyFSPS lib git hash: {pyfsps_lib_hash}")
    hashtag = f"sps_home-{fsps_hash[:6]}_pyfsps-{pyfsps_hash[:6]}_pyfslib-{pyfsps_lib_hash[:6]}"

    sps = fsps.StellarPopulation(zcontinuous=1)
    ilib, slib, dlib = [s.decode("utf-8") for s in sps.libraries]
    print(ilib, slib)
    os.makedirs(f"./figures/{hashtag}", exist_ok=True)
    pdf = PdfPages(f"./figures/{hashtag}/features_{ilib}-{slib}-{dlib}.pdf")

    # Basic spectrum
    sps.params["sfh"] = 4
    sps.params["tau"] = 5.0
    sps.params["logzsol"] = 0.0
    sps.params["dust_type"] = 4  # kriek and Conroy
    sps.params["imf_type"] = 2  # kroupa
    sps.params["imf3"] = 2.3
    fig, ax, spec = makefig(sps)
    fig, ax = prettify(fig, ax, label=r"$\tau=5$, Age$=13.7$,"+"\n" + r"$\log Z/Z_\odot=0.0$")
    pdf.savefig(fig)
    pl.close(fig)

    # change IMF
    sps.params["imf3"] = 2.5
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=f"IMF slope")
    pdf.savefig(fig)

    # Attenuate
    sps.params["add_dust_emission"] = False
    sps.params["dust2"] = 0.2
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=f"Dust Attenuation")
    pdf.savefig(fig)
    pl.close(fig)

    # Dust emission
    sps.params["add_dust_emission"] = True
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=f"Dust Emission")
    pdf.savefig(fig)
    pl.close(fig)

    # Dust temperature
    sps.params["duste_umin"] = 10
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=f"Dust SED\n({dlib})")
    pdf.savefig(fig)
    pl.close(fig)

    # AGN emission
    sps.params["fagn"] = 0.3
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=f"AGN dust\n(Nenkova)")
    pdf.savefig(fig)
    pl.close(fig)

    # Nebular emission
    sps.params["add_neb_emission"] = True
    sps.params["gas_logu"] = -3.5
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=f"Neb. emission\n(Byler)")
    pdf.savefig(fig)
    pl.close(fig)

    # change logu
    sps.params["gas_logu"] = -1.5
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=r"Change U$_{neb}$")
    pdf.savefig(fig)
    pl.close(fig)

    # change logz
    sps.params["logzsol"] = -0.5
    sps.params["gas_logz"] = -0.5
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=r"$\log Z/Z_\odot=-0.5$")
    pdf.savefig(fig)
    pl.close(fig)

    # IGM absorption
    sps.params["zred"] = 6.0
    sps.params["add_igm_absorption"] = True
    fig, ax, spec = makefig(sps, oldspec=spec)
    fig, ax = prettify(fig, ax, label=f"IGM attenuation\n" + r"(Madau, $z=6$)")
    pdf.savefig(fig)
    pl.close(fig)

    pdf.close()
