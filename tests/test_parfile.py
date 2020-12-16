#!/usr/bin/env python3

# Copyright (C) 2020 Gabriele Bozzola
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, see <https://www.gnu.org/licenses/>.

import os

# import pytest

from precactus import parfile


def test_write_one_parfile_from_template():

    outfile = "test_wo.par"

    template = """
################################################################################
# Grid structure
################################################################################

dx = $dx"""
    sub_dict = {"dx": 5}

    parfile.write_one_parfile_from_template(template, sub_dict, outfile)

    with open(outfile, "r") as file_:
        par_file_str = file_.read()

    # We check that the content is what we expect:
    # Variable written by PreCactus
    assert "# dx = 5" in par_file_str
    # Variable substituted in the template
    assert "dx = 5" in par_file_str

    os.remove(outfile)
