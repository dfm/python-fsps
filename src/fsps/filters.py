# -*- coding: utf-8 -*-
"""
Tools for working with the FSPS filter set.

This module uses filter information shipped with FSPS itself in
``$SPS_HOME/data``.
"""

__all__ = ["find_filter", "FILTERS", "get_filter", "list_filters"]

import os
import numpy as np

# Cache for $SPS_HOME/data/magsun.dat parsed by numpy
MSUN_TABLE = None

# Cache for $SPS_HOME/data/filter_lambda_eff.dat parsed by numpy
LAMBDA_EFF_TABLE = None

# Cache for bandpass transmission listings: a dictionary keyed by bandpass
# name with values of wavelength, transmission tuples.
TRANS_CACHE = None


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
        if self.index < MSUN_TABLE.shape[0]:
            assert MSUN_TABLE[self.index, 0] == self.index + 1
            return float(MSUN_TABLE[self.index, 1])
        else:
            return np.nan

    @property
    def msun_vega(self):
        """Solar absolute magnitude in Filter, VEGAMAG zeropoint."""
        if MSUN_TABLE is None:
            self._load_msun_table()
        if self.index < MSUN_TABLE.shape[0]:
            assert MSUN_TABLE[self.index, 0] == self.index + 1
            return float(MSUN_TABLE[self.index, 2])
        else:
            return np.nan

    @property
    def lambda_eff(self):
        """Effective wavelength of Filter, in Angstroms."""
        if LAMBDA_EFF_TABLE is None:
            self._load_lambda_eff_table()
        if self.index < LAMBDA_EFF_TABLE.shape[0]:
            assert LAMBDA_EFF_TABLE[self.index, 0] == self.index + 1
            return float(LAMBDA_EFF_TABLE[self.index, 1])
        else:
            return np.nan

    @property
    def transmission(self):
        """Returns filter transmission: a tuple of wavelength (Angstroms) and
        an un-normalized transmission arrays.
        """
        if TRANS_CACHE is None:
            # Load the cache for all filters.
            self._load_transmission_cache()
        try:
            return TRANS_CACHE[self.name]
        except KeyError as e:
            e.args += "Could not find transmission data " "for {0}".format(self.name)
            raise

    def _load_msun_table(self):
        global MSUN_TABLE
        MSUN_TABLE = np.loadtxt(os.path.expandvars("$SPS_HOME/data/magsun.dat"))

    def _load_lambda_eff_table(self):
        global LAMBDA_EFF_TABLE
        LAMBDA_EFF_TABLE = np.loadtxt(
            os.path.expandvars("$SPS_HOME/data/filter_lambda_eff.dat")
        )

    def _load_transmission_cache(self):
        """Parse the allfilters.dat file into the TRANS_CACHE."""
        global TRANS_CACHE
        path = os.path.expandvars("$SPS_HOME/data/allfilters.dat")
        names = list_filters()
        TRANS_CACHE = {}
        filter_index = -1
        lambdas, trans = [], []
        with open(path) as f:
            for line in f:
                line.strip()
                if line[0].startswith("#"):
                    # Close out filter
                    if filter_index > -1:
                        TRANS_CACHE[names[filter_index]] = (
                            np.array(lambdas),
                            np.array(trans),
                        )
                    # Start new filter
                    filter_index += 1
                    lambdas, trans = [], []
                else:
                    try:
                        l, t = line.split()
                        lambdas.append(float(l))
                        trans.append(float(t))
                    except (ValueError):
                        pass
        # Close out the last filter
        TRANS_CACHE[names[filter_index]] = (np.array(lambdas), np.array(trans))


def _load_filter_dict():
    """
    Load the filter list, creating a dictionary of :class:`Filter` instances.
    """
    # Load filter table from FSPS
    filter_list_path = os.path.expandvars(
        os.path.join("$SPS_HOME", "data", "FILTER_LIST")
    )
    filters = {}
    with open(filter_list_path) as f:
        for line in f:
            columns = line.strip().split()
            if len(columns) > 1:
                fsps_id, key = columns[:2]
            else:
                continue
            comment = " ".join(columns[2:])
            filters[key.lower()] = Filter(int(fsps_id), key, comment)

    return filters


FILTERS = _load_filter_dict()


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


def get_filter(name):
    """Returns the :class:`fsps.filters.Filter` instance associated with the
    filter name.

    :param name:
        Name of the filter, as found with :func:`find_filter`.
    """
    try:
        return FILTERS[name.lower()]
    except KeyError as e:
        e.args += (
            "Filter {0} does not exist. "
            "Try using fsps.find_filter('{0}').".format(name),
        )
        raise


def list_filters():
    """Returns a list of all FSPS filter names.

    Filters are sorted by their FSPS index.
    """
    lst = [(name, f.index) for name, f in FILTERS.items()]
    lst.sort(key=lambda x: x[1])
    return [l[0] for l in lst]
