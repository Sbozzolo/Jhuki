#!/usr/bin/env python3

# Copyright (C) 2021-2022 Gabriele Bozzola
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

import jhuki.simfactory as jsf


def test_simfactory_option():

    assert jsf.simfactory_option("@BOB@", 10) == 10
    # Check int -> float
    assert isinstance(jsf.simfactory_option("@BOB@", 10), float) is True

    assert jsf.simfactory_option("@BOB@", "cat") == "cat"

    assert jsf.simfactory_option("10.4", 20) == 10.4

    # Check wrong type
    with pytest.raises(TypeError):
        assert jsf.simfactory_option(10.4, 20) == 10.4


def test_write_parfile():

    test_var = "HEY"  # noqa: F841

    with pytest.raises(RuntimeError):
        assert jsf.write_parfile("$test_var", "wrong_filename")

    jsf.write_parfile("$test_var", "file.rpar")

    out_file = "file.par"
    with open(out_file, "r") as file_:
        par_file_str = file_.read()

    # We check that the content is what we expect:
    # Variable written by Jhuki
    assert "# test_var = HEY" in par_file_str

    os.remove(out_file)
