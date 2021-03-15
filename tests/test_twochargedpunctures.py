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

from pytest import approx

import jhuki.twochargedpunctures as jtcp


def test_TwoChargedPunctures():

    tp = jtcp.TwoChargedPunctures(
        0.5, 0.5, 12, charge_plus=0.2, charge_minus=-0.2
    )

    assert tp.charge_plus == approx(0.2)
    assert tp.charge_minus == approx(-0.2)

    expected_str = """\
ADMBase::initial_data = "twochargedpunctures"
ADMBase::initial_lapse = "twochargedpunctures-averaged"
ADMBase::initial_shift = "zero"
ADMBase::initial_dtlapse = "zero"
ADMBase::initial_dtshift = "zero"

TwoChargedPunctures::give_bare_mass = "no"
TwoChargedPunctures::par_b = 6.0
TwoChargedPunctures::target_m_plus = 0.5
TwoChargedPunctures::target_m_minus = 0.5
TwoChargedPunctures::par_m_plus = 0.5
TwoChargedPunctures::par_m_minus = 0.5
TwoChargedPunctures::par_P_plus[0] = 0.0
TwoChargedPunctures::par_P_plus[1] = 0.0
TwoChargedPunctures::par_P_plus[2] = 0.0
TwoChargedPunctures::par_P_minus[0] = 0.0
TwoChargedPunctures::par_P_minus[1] = 0.0
TwoChargedPunctures::par_P_minus[2] = 0.0
TwoChargedPunctures::par_S_plus[0] = 0.0
TwoChargedPunctures::par_S_plus[1] = 0.0
TwoChargedPunctures::par_S_plus[2] = 0.0
TwoChargedPunctures::par_S_minus[0] = 0.0
TwoChargedPunctures::par_S_minus[1] = 0.0
TwoChargedPunctures::par_S_minus[2] = 0.0
TwoChargedPunctures::center_offset[0] = 0.0
TwoChargedPunctures::center_offset[1] = 0.0
TwoChargedPunctures::center_offset[2] = 0.0
TwoChargedPunctures::par_q_plus = 0.2
TwoChargedPunctures::par_q_minus = -0.2\
"""

    assert tp.parfile_code == expected_str
