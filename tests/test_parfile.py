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

import pytest

from precactus import parfile


@pytest.fixture(scope="module")
def template():
    template = """
################################################################################
# Grid structure
################################################################################

dx = $dx
dy = $dy"""
    return template


def test_write_one_parfile_from_template(template):

    out_file = "test_wo.par"

    sub_dict = {"dx": 5, "dy": 6}

    parfile.write_one_parfile_from_template(template, sub_dict, out_file)

    with open(out_file, "r") as file_:
        par_file_str = file_.read()

    # We check that the content is what we expect:
    # Variable written by PreCactus
    assert "# dx = 5" in par_file_str
    # Variable substituted in the template
    assert "dx = 5" in par_file_str

    os.remove(out_file)


def test_write_many_parfiles_from_template(template):

    out_file_prefix = "temp_"

    # sub_dict_list contains lists with different lengths
    with pytest.raises(ValueError):
        parfile.write_many_parfiles_from_template(
            template, {"dx": [1, 2], "dy": [3]}, out_file_prefix
        )

    sub_dict = {"dx": [5, 6], "dy": 1}

    parfile.write_many_parfiles_from_template(
        template, sub_dict, out_file_prefix
    )

    # We expect two outfiles:
    for out_file, expected_dx in zip(
        [out_file_prefix + "0.par", out_file_prefix + "1.par"], [5, 6]
    ):
        with open(out_file, "r") as file_:
            par_file_str = file_.read()

        # We check that the content is what we expect:
        assert f"dx = {expected_dx}" in par_file_str
        assert "dy = 1" in par_file_str

        os.remove(out_file)
