#!/usr/bin/env python
# encoding: utf-8
"""
Print a table of filters implemented by FSPS.

The table is formatted in restructured text for use with the python-fsps's
Sphinx documentation.
"""

from __future__ import print_function, absolute_import

import re
import fsps

# Pattern via https://gist.github.com/uogbuji/705383
URL_P = re.compile(ur'(?i)\b((?:https?://|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}/)(?:[^\s()<>]+|\(([^\s()<>]+|(\([^\s()<>]+\)))*\))+(?:\(([^\s()<>]+|(\([^\s()<>]+\)))*\)|[^\s`!()\[\]{};:\'".,<>?\xab\xbb\u201c\u201d\u2018\u2019]))')  # NOQA


def main():
    colnames = ("ID", "Name", "M_sun Vega", "M_sun AB", "lambda_eff (Å)",
                "Description")
    filters = [fsps.get_filter(name) for name in fsps.list_filters()]
    filter_list = make_filter_list(filters)
    txt = make_table(filter_list, colnames)
    print(txt)


def make_filter_list(filters):
    """Transform filters into list of table rows."""
    filter_list = []
    filter_ids = []
    for f in filters:
        filter_ids.append(f.index)
        fullname = URL_P.sub(r'`<\1>`_', f.fullname)
        filter_list.append((str(f.index + 1),
                            f.name,
                            "{0:.2f}".format(f.msun_vega),
                            "{0:.2f}".format(f.msun_ab),
                            "{0:.1f}".format(f.lambda_eff),
                            fullname))
    sortf = lambda item: int(item[0])
    filter_list.sort(key=sortf)
    return filter_list


def make_table(data, col_names):
    """Code for this RST-formatted table generator comes from
    http://stackoverflow.com/a/11350643
    """
    n_cols = len(data[0])
    assert n_cols == len(col_names)
    col_sizes = [max(len(r[i]) for r in data) for i in range(n_cols)]
    for i, cname in enumerate(col_names):
        if col_sizes[i] < len(cname):
            col_sizes[i] = len(cname)
    formatter = ' '.join('{:<%d}' % c for c in col_sizes)
    rows = '\n'.join([formatter.format(*row) for row in data])
    header = formatter.format(*col_names)
    divider = formatter.format(*['=' * c for c in col_sizes])
    output = '\n'.join((divider, header, divider, rows, divider))
    return output


if __name__ == '__main__':
    main()
