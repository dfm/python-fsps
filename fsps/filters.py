#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools for working with the FSPS filter set.

This module uses filter information shipped with FSPS itself in
``$SPS_HOME/data``. Filter names (dictionary keys) are unique to
``python-fsps`` and are defined in ``fsps/data/filter_keys.txt``.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["find_filter", "FILTERS"]

import os
import numpy as np
from pkg_resources import resource_stream, resource_exists

# Cache for $SPS_HOME/data/magsun.dat parsed by numpy
MSUN_TABLE = None

# Cache for $SPS_HOME/data/filter_lambda_eff.dat parsed by numpy
LAMBDA_EFF_TABLE = None


class Filter(object):

    def __init__(self, index, name, fullname):
        self.index = index - 1
        self.name = name.lower()
        self.fullname = fullname

    def __str__(self):
        return "<Filter({0})>".format(self.name)

    def __repr__(self):
        return "<Filter({0})>".format(self.name)

    @property
    def msun_ab(self):
        """Solar absolute magnitude in Filter, AB zeropoint."""
        # if self._msun_ab is None:
        if MSUN_TABLE is None:
            self._load_msun_table()
        return float(MSUN_TABLE[self.index, 1])

    @property
    def msun_vega(self):
        """Solar absolute magnitude in Filter, VEGAMAG zeropoint."""
        if MSUN_TABLE is None:
            self._load_msun_table()
        return float(MSUN_TABLE[self.index, 2])

    @property
    def lambda_eff(self):
        """Effective wavelength of Filter, in Angstroms."""
        if LAMBDA_EFF_TABLE is None:
            self._load_lambda_eff_table()
        return float(LAMBDA_EFF_TABLE[self.index, 1])

    def _load_msun_table(self):
        global MSUN_TABLE
        MSUN_TABLE = np.loadtxt(
            os.path.expandvars("$SPS_HOME/data/magsun.dat"))

    def _load_lambda_eff_table(self):
        global LAMBDA_EFF_TABLE
        LAMBDA_EFF_TABLE = np.loadtxt(
            os.path.expandvars("$SPS_HOME/data/filter_lambda_eff.dat"))


def load_filter_dict():
    """
    Load the filter list, creating a dictionary of :class:`Filter` instances.
    """
    # Load filter table from FSPS
    filter_list_path = os.path.expandvars(
        os.path.join("$SPS_HOME", "data", "FILTER_LIST"))
    fsps_ids, comments = [], []
    with open(filter_list_path) as f:
        for line in f:
            fsps_id, comment = line.split('\t')
            fsps_ids.append(int(fsps_id))
            comments.append(comment)

    # Load our own table of filter keys
    our_ids, keys = [], []
    keys_path = "data/filter_keys.txt"
    assert resource_exists(__name__, keys_path)
    for line in resource_stream(__name__, keys_path):
        _id, key = line.split()
        our_ids.append(int(_id))
        keys.append(key)

    # Merge tables
    assert len(fsps_ids) == len(our_ids)
    filters = {}
    for fsps_id, filter_key, comment in zip(fsps_ids, keys, comments):
        filters[filter_key.lower()] = Filter(fsps_id, filter_key, comment)
    return filters


FILTERS = load_filter_dict()


def find_filter(band):
    """
    Find the FSPS name for a filter.

    Usage:

    ::

        >>> import fsps
        >>> fsps.find_filter("F555W")
        ['wfpc2_f555w', 'wfc_acs_f555w']

    :param band:
        Something like the name of the band.

    """
    b = band.lower()
    possible = []
    for k in FILTERS.keys():
        if b in k:
            possible.append(k)
    return possible
