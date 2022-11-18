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
from pytest import approx

import jhuki.twopunctures as jtp


def test_TwoPunctures():

    # Negative masses
    with pytest.raises(ValueError):
        # Plus
        jtp.TwoPunctures(-0.5, 0.5, 12)
    with pytest.raises(ValueError):
        # Minus
        jtp.TwoPunctures(0.5, -0.5, 12)

    # Both zero masses
    with pytest.raises(ValueError):
        jtp.TwoPunctures(0.0, -0.0, 12)

    # Mass-ratio 1, everything else default
    tp = jtp.TwoPunctures(0.5, 0.5, 12)

    assert tp.par_b == approx(6)
    assert tp.coord_x_plus == approx(6)
    assert tp.coord_x_minus == approx(-6)

    assert tp.momenta_plus == (0, 0, 0)
    assert tp.momenta_minus == (0, 0, 0)
    assert tp.chi_plus == (0, 0, 0)
    assert tp.chi_minus == (0, 0, 0)
    assert tp.S_plus == (0, 0, 0)
    assert tp.S_minus == (0, 0, 0)
    assert tp.center_offset == (0, 0, 0)

    expected_str = """\
ADMBase::initial_data = "twopunctures"
ADMBase::initial_lapse = "twopunctures-averaged"
ADMBase::initial_shift = "zero"
ADMBase::initial_dtlapse = "zero"
ADMBase::initial_dtshift = "zero"

TwoPunctures::give_bare_mass = "no"
TwoPunctures::par_b = 6.0
TwoPunctures::target_m_plus = 0.5
TwoPunctures::target_m_minus = 0.5
TwoPunctures::par_m_plus = 0.5
TwoPunctures::par_m_minus = 0.5
TwoPunctures::par_P_plus[0] = 0.0
TwoPunctures::par_P_plus[1] = 0.0
TwoPunctures::par_P_plus[2] = 0.0
TwoPunctures::par_P_minus[0] = 0.0
TwoPunctures::par_P_minus[1] = 0.0
TwoPunctures::par_P_minus[2] = 0.0
TwoPunctures::par_S_plus[0] = 0.0
TwoPunctures::par_S_plus[1] = 0.0
TwoPunctures::par_S_plus[2] = 0.0
TwoPunctures::par_S_minus[0] = 0.0
TwoPunctures::par_S_minus[1] = 0.0
TwoPunctures::par_S_minus[2] = 0.0
TwoPunctures::center_offset[0] = 0.0
TwoPunctures::center_offset[1] = 0.0
TwoPunctures::center_offset[2] = 0.0\
"""

    assert tp.parfile_code == expected_str
    assert str(tp) == tp.parfile_code

    # Mass-ratio 2, everything else bizarre
    tp2 = jtp.TwoPunctures(
        1,
        2,
        18,
        (0.1, 0.2, 0.3),
        (0.4, 0.5, 0.6),
        (-0.1, -0.2, -0.3),
        (-0.4, -0.5, -0.6),
        swap_xz=True,
        give_bare_mass=True,
        initial_alpha="psi^n",
    )

    assert tp2.par_b == approx(9)
    assert tp2.coord_x_plus == approx(12)
    assert tp2.coord_x_minus == approx(-6)

    assert tp2.momenta_plus == (0.1, 0.2, 0.3)
    assert tp2.momenta_minus == (0.4, 0.5, 0.6)
    assert tp2.chi_plus == (-0.1, -0.2, -0.3)
    assert tp2.chi_minus == (-0.4, -0.5, -0.6)
    assert tp2.a_plus == (-0.1, -0.2, -0.3)
    assert tp2.a_minus == (-0.8, -1.0, -1.2)
    assert tp2.S_plus == (-0.1, -0.2, -0.3)
    assert tp2.S_minus == (-1.6, -2.0, -2.4)
    assert tp2.center_offset == (3, 0, 0)

    # Note, we have to swap x z in the vectors!
    expected_str = """\
ADMBase::initial_data = "twopunctures"
ADMBase::initial_lapse = "psi^n"
ADMBase::initial_shift = "zero"
ADMBase::initial_dtlapse = "zero"
ADMBase::initial_dtshift = "zero"

TwoPunctures::give_bare_mass = "yes"
TwoPunctures::par_b = 9.0
TwoPunctures::target_m_plus = 1
TwoPunctures::target_m_minus = 2
TwoPunctures::par_m_plus = 1
TwoPunctures::par_m_minus = 2
TwoPunctures::par_P_plus[0] = 0.3
TwoPunctures::par_P_plus[1] = 0.2
TwoPunctures::par_P_plus[2] = 0.1
TwoPunctures::par_P_minus[0] = 0.6
TwoPunctures::par_P_minus[1] = 0.5
TwoPunctures::par_P_minus[2] = 0.4
TwoPunctures::par_S_plus[0] = -0.3
TwoPunctures::par_S_plus[1] = -0.2
TwoPunctures::par_S_plus[2] = -0.1
TwoPunctures::par_S_minus[0] = -2.4
TwoPunctures::par_S_minus[1] = -2.0
TwoPunctures::par_S_minus[2] = -1.6
TwoPunctures::center_offset[0] = 0.0
TwoPunctures::center_offset[1] = 0.0
TwoPunctures::center_offset[2] = 3.0
TwoPunctures::swap_xz = yes\
"""

    assert tp2.parfile_code == expected_str

    # One mass zero, everything else default
    tp3 = jtp.TwoPunctures(0.0, 0.5, 12)
    assert tp3.par_b == approx(6)
    assert tp3.coord_x_plus == approx(12)
    assert tp3.coord_x_minus == approx(0)
    assert tp3.center_offset == (6, 0, 0)


def test_prepare_quasicircular_inspiral():

    # Test mass ratio < 1
    tp3 = jtp.prepare_quasicircular_inspiral(0.2, 12)

    assert tp3.mass_plus == 5 / 6
    assert tp3.mass_minus == 1 / 6

    assert tp3.momenta_plus == approx(
        (-0.000169968781552016, 0.0474161839456146, 0)
    )
    assert tp3.momenta_minus == approx(
        (0.000169968781552016, -0.0474161839456146, 0)
    )

    # Test mass ratio > 1
    tp4 = jtp.prepare_quasicircular_inspiral(5, 12)

    assert tp4.mass_plus == 5 / 6
    assert tp4.mass_minus == 1 / 6

    assert tp4.momenta_plus == approx(
        (-0.000169968781552016, 0.0474161839456146, 0)
    )
    assert tp4.momenta_minus == approx(
        (0.000169968781552016, -0.0474161839456146, 0)
    )

    # Test total_bare_mass = 10
    total_bare_mass = 10
    tp5 = jtp.prepare_quasicircular_inspiral(
        0.2, 12, total_bare_mass=total_bare_mass
    )

    assert tp5.mass_plus == 5 / 6 * total_bare_mass
    assert tp5.mass_minus == 1 / 6 * total_bare_mass
    assert tp5.momenta_plus == approx(
        (
            -0.000169968781552016 * total_bare_mass,
            0.0474161839456146 * total_bare_mass,
            0,
        )
    )
