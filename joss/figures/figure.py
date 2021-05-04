#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pl
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
        ax.set_title(label, fontsize=16)
    return fig, ax


pl.rc("axes", grid=False)
pl.rc("xtick", direction="in")
pl.rc("ytick", direction="in")
pl.rc("xtick", top=True)
pl.rc("ytick", right=True)

sps = fsps.StellarPopulation(zcontinuous=1)
ilib, slib, dlib = sps.libraries

# Basic spectrum
sps.params["sfh"] = 4
sps.params["tau"] = 5.0
sps.params["logzsol"] = 0.0
sps.params["dust_type"] = 4  # kriek and Conroy
sps.params["imf_type"] = 2  # kroupa
sps.params["imf3"] = 2.3
fig, ax, spec = makefig(sps)
fig, ax = prettify(
    fig, ax, label=r"$\tau=5$; Age$=13.7$; $\log Z/Z_\odot=0.0$"
)

pl.savefig("figure.png", dpi=200, bbox_inches="tight")
