#!/usr/bin/env python3

# Copyright (C) 2021 Gabriele Bozzola
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

from jhuki import dissipation as dis
from jhuki import grid as gr


@pytest.fixture(scope="module")
def ds():
    ds = dis.Dissipation(
        {0: 0.1, 1: 0.1, 2: 0.2},
        order=5,
        variables=["Thorn1::var1", "Thorn2::var2"],
    )
    return ds


def test_parfile_code(ds):

    expected_str = """\
Dissipation::order = 5
Dissipation::vars  = "Thorn1::var1
Thorn2::var2"
Dissipation::epsdis_for_level[0] = 0.1
Dissipation::epsdis_for_level[1] = 0.1
Dissipation::epsdis_for_level[2] = 0.2\
"""

    assert ds.parfile_code == expected_str


def test_create_dissipation_from_grid():

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

    # dtfact
    ds = dis.create_dissipation_from_grid(grid, 0.1)

    assert ds.epsdis_per_level == {
        0: 0.025,
        1: 0.05,
        2: 0.1,
        3: 0.1,
        4: 0.1,
        5: 0.1,
        6: 0.1,
    }

    # constant
    ds = dis.create_dissipation_from_grid(
        grid, 0.1, prescription=dis.DissPrescription.const
    )

    assert ds.epsdis_per_level == {
        0: 0.1,
        1: 0.1,
        2: 0.1,
        3: 0.1,
        4: 0.1,
        5: 0.1,
        6: 0.1,
    }

    # continuous
    ds = dis.create_dissipation_from_grid(
        grid, 0.1, prescription=dis.DissPrescription.continuous
    )

    assert ds.epsdis_per_level == {
        0: 0.1 / 2 ** 30,
        1: 0.1 / 2 ** 25,
        2: 0.1 / 2 ** 20,
        3: 0.1 / 2 ** 15,
        4: 0.1 / 2 ** 10,
        5: 0.1 / 2 ** 5,
        6: 0.1,
    }

    # invalid
    with pytest.raises(ValueError):
        dis.create_dissipation_from_grid(grid, 0.1, prescription="bob")


def test_add_Lean(ds):
    ds_with_Lean = dis.add_Lean(ds)

    expected_vars = ds.variables + [
        "LeanBSSNMoL::conf_fac",
        "LeanBSSNMoL::hmetric",
        "LeanBSSNMoL::hcurv",
        "LeanBSSNMoL::trk",
        "LeanBSSNMoL::gammat",
    ]

    assert ds_with_Lean.variables == expected_vars


def test_add_gauge(ds):
    ds_with_gauge = dis.add_gauge(ds)

    expected_vars = ds.variables + ["ADMBase::lapse", "ADMBase::shift"]

    assert ds_with_gauge.variables == expected_vars


def test_add_Proca(ds):
    ds_with_Proca = dis.add_Proca(ds)

    expected_vars = ds.variables + [
        "ProcaBase::Ei",
        "ProcaBase::Ai",
        "ProcaBase::Aphi",
        "ProcaBase::Zeta",
    ]

    assert ds_with_Proca.variables == expected_vars
