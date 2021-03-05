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

import jhuki.externals.nrpypn as nrpypn

# Check that we can use eval__P_t__and__P_r


def test_compute_quasicircular_momenta():
    assert nrpypn.compute_quasicircular_momenta(
        1, 12, (0.1, 0.2, 0.3), (-0.1, -0.2, -0.3)
    ) == approx((0.0850933618449702, 0.000539545145231776))
