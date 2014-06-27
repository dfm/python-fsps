#!/usr/bin/env python
# encoding: utf-8
"""
Demonstration of plotting FSPS's un-normalized filter transmission tables
(i.e., contents of $SPS_HOME/data/all_filters.dat).
"""

import fsps

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec


names = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i',
         '2mass_J', '2mass_H', '2mass_Ks']
shortnames = ['u', 'g', 'r', 'i', 'J', 'H', 'Ks']
colors = ['violet', 'dodgerblue', 'maroon', 'black', 'c', 'm', 'y']
filters = [fsps.get_filter(n) for n in names]

fig = Figure(figsize=(3.5, 3.5))
canvas = FigureCanvas(fig)
gs = gridspec.GridSpec(
    1, 1, left=0.17, right=0.95, bottom=0.15, top=0.95,
    wspace=None, hspace=None, width_ratios=None, height_ratios=None)
ax = fig.add_subplot(gs[0])

for name, fltr, c, shortname in zip(names, filters, colors, shortnames):
    lmbd, trans = fltr.transmission
    lambda_eff = fltr.lambda_eff / 10000.  # µm
    ax.plot(lmbd / 10000., trans, ls='-', lw=1., c=c)
    ax.annotate(shortname, (lambda_eff, trans.max()),
                textcoords='offset points',
                xytext=(0., 5.),
                ha='center', va='bottom')

ax.set_ylim(0., 1.1)
ax.set_xlabel(r"$\lambda~\mu\mathrm{m}$")
ax.set_ylabel(r"$T$")
gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
canvas.print_figure("sdss_2mass_transmission_demo.pdf", format="pdf")
