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

from math import sqrt

import pytest

from jhuki import grid as gr
from jhuki.twochargedpunctures import TwoChargedPunctures
from jhuki.twopunctures import TwoPunctures


@pytest.fixture(scope="module")
def refinement_center():
    mr = gr.RefinementCenter(
        refinement_radii=[1, 2, 3],
        dx_fine=0.25,
        cfl_fine=1.0,
        center_num=2,
        position=(1, 0, -1),
        num_levels_with_dt_coarse=2,
    )
    return mr


@pytest.fixture(scope="module")
def refinement_center2():
    mr = gr.RefinementCenter(
        refinement_radii=[0.1, 0.2, 0.3, 4, 5],
        dx_fine=0.0625,
        cfl_fine=1,
        center_num=1,
        position=(2, 5, 3),
        num_levels_with_dt_coarse=2,
    )
    return mr


@pytest.fixture(scope="module")
def a_grid(refinement_center, refinement_center2):

    grid = gr.Grid([refinement_center, refinement_center2], outer_boundary=10)

    return grid


def test_refinement_center_init(refinement_center):

    # Refinement_radii not a list
    with pytest.raises(TypeError):
        gr.RefinementCenter("HI", 1, 0.2)

    # Refinement_radii is empty
    with pytest.raises(ValueError):
        gr.RefinementCenter([], 1, 0.2)

    # Invalid num_levels_with_dt_coarse
    with pytest.raises(ValueError):
        gr.RefinementCenter([1, 2], 1, 0.2, num_levels_with_dt_coarse=0)


def test_dx(refinement_center):
    expected_dx = (2.0, 1.0, 0.5, 0.25)

    assert refinement_center.dx == expected_dx


def test_dt(refinement_center):
    expected_dt = (1.0, 1.0, 0.5, 0.25)

    assert refinement_center.dt == expected_dt


def test_rl_synced_every(refinement_center, refinement_center2, a_grid):

    # refinement_center has 3 refinement levels
    # and num_levels_with_dt_coarse=2, so the system is
    # OB - RL0 - RL1 - RL2 - (0,0,0)
    #   dt    dt    dt/2  dt/4
    #
    # So, we have to evolve twice 4 timesteps the innermost level,
    #
    assert refinement_center.rl_synced_every == 4

    # refinement_center2 has 5 levels instead
    assert refinement_center2.rl_synced_every == 16

    # The grid should consider the one with the most levels
    assert a_grid.rl_synced_every == 16


def test_cfl(refinement_center):
    expected_cfl = (0.5, 1.0, 1.0, 1.0)

    assert refinement_center.cfl == expected_cfl


def test_hash(refinement_center):

    hash_refinement_radii = hash(refinement_center.refinement_radii)
    hash_center_num = hash(refinement_center.center_num)
    hash_position = hash(refinement_center.position)
    hash_dx_fine = hash(refinement_center.dx_fine)
    hash_cfl_fine = hash(refinement_center.cfl_fine)
    hash_num_levels_with_dt_coarse = hash(
        refinement_center.num_levels_with_dt_coarse
    )

    expected_hash = (
        hash_refinement_radii
        ^ hash_center_num
        ^ hash_position
        ^ hash_dx_fine
        ^ hash_cfl_fine
        ^ hash_num_levels_with_dt_coarse
    )

    assert expected_hash == hash(refinement_center)


def test_parfile_code(refinement_center):

    expected_str = """\
CarpetRegrid2::num_levels_2 = 4
CarpetRegrid2::position_x_2 = 1
CarpetRegrid2::position_y_2 = 0
CarpetRegrid2::position_z_2 = -1
CarpetRegrid2::radius_2[1] = 3
CarpetRegrid2::radius_2[2] = 2
CarpetRegrid2::radius_2[3] = 1\
"""

    assert refinement_center.parfile_code == expected_str


def test_summary(refinement_center):

    expected_str = """\
Refinement center 2 (x =  1.00, y =  0.00,  z =  0.00)
Ref level  0: Radius range = (    5.000,     3.000), dx = 2.0000e+00, dt = 1.0000e+00, cfl = 0.500
Ref level  1: Radius range = (    3.000,     2.000), dx = 1.0000e+00, dt = 1.0000e+00, cfl = 1.000
Ref level  2: Radius range = (    2.000,     1.000), dx = 5.0000e-01, dt = 5.0000e-01, cfl = 1.000
Ref level  3: Radius range = (    1.000,     0.000), dx = 2.5000e-01, dt = 2.5000e-01, cfl = 1.000\
"""
    assert expected_str == refinement_center.get_summary(outer_boundary=5)

    expected_str9 = """\
Refinement center 2 (x =  1.00, y =  0.00,  z =  0.00)
Ref level  0: Radius range = ( 9999.999,     3.000), dx = 2.0000e+00, dt = 1.0000e+00, cfl = 0.500
Ref level  1: Radius range = (    3.000,     2.000), dx = 1.0000e+00, dt = 1.0000e+00, cfl = 1.000
Ref level  2: Radius range = (    2.000,     1.000), dx = 5.0000e-01, dt = 5.0000e-01, cfl = 1.000
Ref level  3: Radius range = (    1.000,     0.000), dx = 2.5000e-01, dt = 2.5000e-01, cfl = 1.000\
"""
    assert expected_str9 == refinement_center.get_summary()


def test_set_dt_max(refinement_center):

    # Too small dt_max
    with pytest.raises(ValueError):
        gr.set_dt_max(refinement_center, 0.0001)

    with pytest.raises(ValueError):
        gr.set_dt_max(refinement_center, 0.0001)

    assert (
        gr.set_dt_max(refinement_center, 1.20).num_levels_with_dt_coarse == 2
    )
    assert (
        gr.set_dt_max(refinement_center, 0.55).num_levels_with_dt_coarse == 3
    )
    assert (
        gr.set_dt_max(refinement_center, 0.33).num_levels_with_dt_coarse == 4
    )


def test_grid_init(refinement_center, a_grid):

    # Wrong type
    with pytest.raises(TypeError):
        gr.Grid([1], 10)

    # Same center_num
    with pytest.raises(ValueError):
        gr.Grid(
            (
                refinement_center,
                gr.RefinementCenter(
                    refinement_radii=[1],
                    dx_fine=0.25,
                    cfl_fine=1,
                    center_num=2,
                ),
            ),
            outer_boundary=10,
        )

    # Different dxs
    with pytest.raises(ValueError):
        gr.Grid(
            (
                refinement_center,
                gr.RefinementCenter(
                    refinement_radii=[1], dx_fine=0.25, cfl_fine=1
                ),
            ),
            outer_boundary=10,
        )

    # Different dt_coarse
    with pytest.raises(ValueError):
        gr.Grid(
            (
                refinement_center,
                gr.RefinementCenter(
                    refinement_radii=[1, 2, 3],
                    dx_fine=0.25,
                    cfl_fine=10,
                    center_num=1,
                    num_levels_with_dt_coarse=2,
                ),
            ),
            outer_boundary=10,
        )

    # Different num_dt_coarse
    with pytest.raises(ValueError):
        gr.Grid(
            (
                refinement_center,
                gr.RefinementCenter(
                    refinement_radii=[1, 2, 3],
                    dx_fine=0.25,
                    cfl_fine=0.5,
                    center_num=1,
                    num_levels_with_dt_coarse=1,
                ),
            ),
            outer_boundary=10,
        )

    # Compatible refinement centers but with different
    # num_levels_with_dt_coarse.
    griddo = gr.Grid(
        (
            gr.RefinementCenter(
                refinement_radii=[1, 2, 3],
                dx_fine=0.25,
                cfl_fine=1.0,
                center_num=1,
                position=(1, 0, -1),
                num_levels_with_dt_coarse=4,
            ),
            gr.RefinementCenter(
                refinement_radii=[1, 2, 3, 4, 5],
                dx_fine=0.0625,
                cfl_fine=1.0,
                center_num=2,
                position=(0, 0, 0),
                num_levels_with_dt_coarse=4,
            ),
        ),
        outer_boundary=10,
    )

    assert griddo.time_refinement_factors == '"[1,1,1,1,2,4]"'

    # Test _coarse
    assert a_grid.dt_coarse == 1
    assert a_grid.dx_coarse == 2
    assert a_grid.cfl_coarse == 0.5


def test_time_refinement_factors(a_grid):

    expected_time_refinement = '"[1,1,2,4,8,16]"'

    assert expected_time_refinement == a_grid.time_refinement_factors


def test_reflection_axis(refinement_center, refinement_center2):

    with pytest.raises(ValueError):
        gr.Grid(
            [refinement_center, refinement_center2],
            outer_boundary=10,
            reflection_axis="bob",
        )


def test_grid_parfile_code(a_grid, refinement_center, refinement_center2):

    expected_str = f"""\
CartGrid3D::type = "coordbase"
Carpet::domain_from_coordbase = "yes"

Driver::ghost_size = 3
CoordBase::boundary_size_x_lower = 3
CoordBase::boundary_size_y_lower = 3
CoordBase::boundary_size_z_lower = 3
CoordBase::boundary_size_x_upper = 3
CoordBase::boundary_size_y_upper = 3
CoordBase::boundary_size_z_upper = 3

CoordBase::domainsize = "minmax"
CoordBase::xmax = 10
CoordBase::ymax = 10
CoordBase::zmax = 10
CoordBase::xmin = -10
CoordBase::ymin = -10
CoordBase::zmin = -10
CoordBase::dx = 2.0
CoordBase::dy = 2.0
CoordBase::dz = 2.0

CarpetRegrid2::num_centres = 2
Carpet::max_refinement_levels = 6
Carpet::time_refinement_factors = "[1,1,2,4,8,16]"

{a_grid.refinement_centers[0].parfile_code}
{a_grid.refinement_centers[1].parfile_code}\
"""

    assert a_grid.parfile_code == expected_str

    grid_symmetry_x = gr.Grid(
        [refinement_center, refinement_center2],
        outer_boundary=10,
        reflection_axis="x",
        num_ghost=4,
    )

    expected_str = f"""\
CartGrid3D::type = "coordbase"
Carpet::domain_from_coordbase = "yes"

Driver::ghost_size = 4
CoordBase::boundary_size_x_lower = 4
CoordBase::boundary_size_y_lower = 4
CoordBase::boundary_size_z_lower = 4
CoordBase::boundary_size_x_upper = 4
CoordBase::boundary_size_y_upper = 4
CoordBase::boundary_size_z_upper = 4

CoordBase::domainsize = "minmax"
CoordBase::xmax = 10
CoordBase::ymax = 10
CoordBase::zmax = 10
CoordBase::xmin = 0
CoordBase::ymin = -10
CoordBase::zmin = -10
CoordBase::dx = 2.0
CoordBase::dy = 2.0
CoordBase::dz = 2.0

ReflectionSymmetry::reflection_x     = yes
ReflectionSymmetry::avoid_origin_x   = no
CoordBase::boundary_shiftout_x_lower = 1

CarpetRegrid2::num_centres = 2
Carpet::max_refinement_levels = 6
Carpet::time_refinement_factors = "[1,1,2,4,8,16]"

{a_grid.refinement_centers[0].parfile_code}
{a_grid.refinement_centers[1].parfile_code}\
"""

    assert grid_symmetry_x.parfile_code == expected_str


def test_set_dt_max_grid(a_grid):

    mr1 = gr.RefinementCenter(
        refinement_radii=[1, 2, 3],
        dx_fine=0.25,
        cfl_fine=1.0,
        center_num=2,
        position=(1, 0, -1),
        num_levels_with_dt_coarse=3,
    )

    mr2 = gr.RefinementCenter(
        refinement_radii=[0.1, 0.2, 0.3, 4, 5],
        dx_fine=0.0625,
        cfl_fine=1,
        center_num=1,
        position=(2, 5, 3),
        num_levels_with_dt_coarse=3,
    )

    expected_grid = gr.Grid((mr1, mr2), outer_boundary=10)

    # We haven't defined the __eq__ method, so we check that the parfile codes
    # are the same
    assert (
        expected_grid.parfile_code
        == gr.set_dt_max_grid(a_grid, 0.5).parfile_code
    )

    # Test with tiny_shift and reflection symmetry
    grido = gr.Grid(
        (mr1,),
        outer_boundary=5,
        tiny_shift=True,
        num_ghost=4,
        reflection_axis="xz",
    )

    outer_boundary_plus = 5 + 0.25 / 7

    assert str(outer_boundary_plus) in grido.parfile_code

    assert "ReflectionSymmetry::reflection_z     = yes" in grido.parfile_code

    assert grido.num_ghost == 4

    # Test with one case in which an entire refinement center has to be synced
    # up

    mr1_2 = gr.RefinementCenter(
        refinement_radii=[1, 2, 3],
        dx_fine=0.25,
        cfl_fine=1.0,
        center_num=1,
        position=(1, 0, -1),
        num_levels_with_dt_coarse=2,
    )
    mr2_2 = gr.RefinementCenter(
        refinement_radii=[1, 2, 3, 4, 5],
        dx_fine=0.0625,
        cfl_fine=1.0,
        center_num=2,
        position=(0, 0, 0),
        num_levels_with_dt_coarse=2,
    )

    griddo = gr.set_dt_max_grid(
        gr.Grid((mr1_2, mr2_2), outer_boundary=10), 0.13
    )

    assert griddo.time_refinement_factors == '"[1,1,1,1,1,2]"'


def test_create_twopuncture_grid():

    tp = TwoPunctures(
        mass_plus=0.5,
        mass_minus=0.3,
        coordinate_distance=10,
        chi_plus=(0.1, 0.2, 0.3),
        chi_minus=(0.15, 0.25, 0.35),
    )

    grid_tp = gr.create_twopunctures_grid(
        tp,
        points_on_horizon_radius=40,
        minimum_outer_boundary=100,
        tiny_shift=True,
    )

    # We cover the mass_minus horizon, which is the smallest.

    spin_m = 0.3 * sqrt(0.15 ** 2 + 0.25 ** 2 + 0.35 ** 2)

    ah_m = sqrt(0.3 ** 2 - spin_m ** 2) / 2

    dx_fine = round(ah_m / 40, 4)

    assert grid_tp.tiny_shift is True
    assert grid_tp.refinement_centers[0].dx_fine == dx_fine
    assert grid_tp.refinement_centers[0].position == (tp.coord_x_plus, 0, 0)
    assert grid_tp.refinement_centers[1].position == (tp.coord_x_minus, 0, 0)
    assert grid_tp.refinement_centers[0].num_refinement_radii == 10
    assert grid_tp.outer_boundary == pytest.approx(152.064)

    tcp = TwoChargedPunctures(
        mass_plus=0.5,
        mass_minus=0.3,
        coordinate_distance=10,
        chi_plus=(0.1, 0.2, 0.3),
        chi_minus=(0.15, 0.25, 0.35),
        charge_minus=0.03,
        swap_xz=True,
    )

    grid_tcp = gr.create_twopunctures_grid(
        tcp, points_on_horizon_radius=40, minimum_outer_boundary=100
    )

    assert grid_tcp.refinement_centers[0].position == (0, 0, tcp.coord_x_plus)
    assert grid_tcp.refinement_centers[1].position == (0, 0, tcp.coord_x_minus)

    ah_m = sqrt(0.3 ** 2 - spin_m ** 2 - 0.03 ** 2) / 2
    dx_fine = round(ah_m / 40, 4)
    assert grid_tcp.refinement_centers[0].dx_fine == dx_fine
