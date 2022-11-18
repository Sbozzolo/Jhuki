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

import pytest

from jhuki import grid as gr
from jhuki import nanchecker as nan


@pytest.fixture(scope="module")
def nanc():
    nanc = nan.NaNChecker(
        variables=["Thorn1::var1", "Thorn2::var2"],
        check_every=128,
        action_if_found="abort",
    )
    return nanc


def test_parfile_code(nanc):

    expected_str = """\
NaNChecker::check_every = 128
NaNChecker::action_if_found = abort
NaNChecker::check_vars  = "Thorn1::var1
Thorn2::var2"\
"""

    assert nanc.parfile_code == expected_str


def test_create_nanchecker_from_grid(nanc):

    mr1 = gr.RefinementCenter(
        refinement_radii=[1, 2, 3],
        dx_fine=0.25,
        cfl_fine=1.0,
        center_num=2,
        position=(1, 0, -1),
        num_levels_with_dt_coarse=2,
    )

    mr2 = gr.RefinementCenter(
        refinement_radii=[0.1, 0.2, 0.3, 4, 5],
        dx_fine=0.0625,
        cfl_fine=1,
        center_num=1,
        position=(2, 5, 3),
        num_levels_with_dt_coarse=2,
    )

    grid = gr.Grid([mr1, mr2], outer_boundary=10)

    nanny = nan.create_nanchecker_from_grid(grid, nanc.variables)

    assert nanny.check_every == 16
