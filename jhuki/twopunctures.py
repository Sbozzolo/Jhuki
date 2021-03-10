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

"""The :py:mod:`~.twopunctures` module provides infrastructure to prepare
binary black simulations.

Some functions are provided:
- :py:func:`~.prepare_quasicircular_inspiral`, which returns a :py:class:`~.TwoPunctures`
  with parameters set up to lead to a quasi-circular inspiral. The linear momenta
  required for this are computed with a 3.5PN expansion by NRPyPN.

"""

from functools import lru_cache

from jhuki.externals.nrpypn import compute_quasicircular_momenta

# TODO: This module needs unit testing


def prepare_quasicircular_inspiral(
    mass_ratio, coordinate_distance, chi_plus=(0, 0, 0), chi_minus=(0, 0, 0)
):
    """Return a :py:class:`~.TwoPunctures` that describes a quasi-circular inspiral.
    We always assume that the plus puncture is the most massive one.

    :param mass_ratio: Mass ratio of the binary.
    :type mass_ratio: float
    :param coordinate_distance: Initial coordinate separation.
    :type coordinate_distance: float

    :param chi_plus: Dimensionless spin of the black hole on the positive side of the x
                     (or z) axis along the three directions.
    :type chi_plus: tuple/list with three numbers

    :param chi_minus: Dimensionless spin of the black hole on the negative side of the x
                     (or z) axis along the three directions.
    :type chi_minus: tuple/list with three numbers

    :returns: A :py:class:`~.TwoPunctures` for a quasi-circular inspiral.
    :rtype: :py:class:`~.TwoPunctures`
    """

    if mass_ratio < 1:
        mass_ratio = 1 / mass_ratio

    Pt, Pr = compute_quasicircular_momenta(
        mass_ratio, coordinate_distance, chi_plus, chi_minus
    )

    momenta_plus = (-Pr, Pt, 0)
    momenta_minus = (Pr, -Pt, 0)

    mass_plus = mass_ratio / (1 + mass_ratio)
    mass_minus = 1 / (1 + mass_ratio)

    return TwoPunctures(
        mass_plus,
        mass_minus,
        coordinate_distance,
        momenta_plus=momenta_plus,
        momenta_minus=momenta_minus,
        chi_plus=chi_plus,
        chi_minus=chi_minus,
    )


class TwoPunctures:
    """The :py:class:`~.TwoPunctures` class contains all the information related to
    one or two black holes, their locations, and properties.

    This class can set up parameters for the ``TwoPunctures`` thorn making sure that
    the center of mass is in the center of the grid (unless differently specified).

    This class is immutable.

    :ivar mass_plus: Mass of the black hole on the positive side of the x (or z) axis.
    :vartype mass_plus: float
    :ivar mass_minus: Mass of the black hole on the negative side of the x (or z) axis.
    :vartype mass_minus: float
    :ivar coordinate_separation: Initial coordinate separation.
    :vartype coordinate_separation: float

    :ivar coord_x_plus: Initial coordinate location of the center of the black hole on the
                        positive side of the x (or z) axis.
    :vartype coord_x_plus: float
    :ivar coord_x_minus: Initial coordinate location of the center of the black hole on the
                         negative side of the x (or z) axis.
    :vartype coord_x_minus: float

    :ivar momenta_plus: Array with the linear momenta along the three directions for the
                        black hole on the positive side of the x (or z) axis.
    :vartype momenta_plus: list/tuple
    :ivar momenta_minus: Array with the linear momenta along the three directions for the
                         black hole on the positive side of the x (or z) axis.
    :vartype momenta_minus: list/tuple

    :ivar chi_plus: Array with the dimensionless_spin along the three directions for the
                    black hole on the positive side of the x (or z) axis.
    :vartype chi_plus: list/tuple
    :ivar chi_minus: Array with the dimensionless spin along the three directions for the
                     black hole on the positive side of the x (or z) axis.
    :vartype chi_minus: list/tuple

    :ivar S_plus: Array with the angular momentum along the three directions for the
                  black hole on the positive side of the x (or z) axis.
    :vartype S_plus: list/tuple
    :ivar S_minus: Array with the angular momentum along the three directions for the
                   black hole on the positive side of the x (or z) axis.
    :vartype S_minus: list/tuple

    :ivar center_offset: Move the geometric center of the system by this amount.
    :vartype center_offset: list/tuple
    """

    def __init__(
        self,
        mass_plus,
        mass_minus,
        coordinate_distance,
        momenta_plus=None,
        momenta_minus=None,
        chi_plus=None,
        chi_minus=None,
        center_offset=None,
    ):
        """Constructor.

        :param mass_plus: Mass of the black hole on the positive side of the x (or z) axis.
        :type mass_plus: float
        :param mass_minus: Mass of the black hole on the negative side of the x (or z) axis.
        :type mass_minus: float
        :param coordinate_separation: Initial coordinate separation.
        :type coordinate_separation: float
        :param par_b: Half of the initial separation.
        :type par_b: float

        :param momenta_plus: Array with the linear momenta along the three directions for the
                            black hole on the positive side of the x (or z) axis.
        :type momenta_plus: list/tuple
        :param momenta_minus: Array with the linear momenta along the three directions for the
                             black hole on the positive side of the x (or z) axis.
        :type momenta_minus: list/tuple

        :param chi_plus: Array with the dimensionless_spin along the three directions for the
                        black hole on the positive side of the x (or z) axis.
        :type chi_plus: list/tuple
        :param chi_minus: Array with the dimensionless spin along the three directions for the
                         black hole on the positive side of the x (or z) axis.
        :type chi_minus: list/tuple

        :ivar center_offset: Move the center of mass of the system by this amount. If None, it
                             this is set so that the center of mass is at [0,0,0].
        :vartype center_offset: list/tuple
        """

        # TODO: Support swapping x and z axis
        # TODO: No sanity checks are performed here. We should at least check that P < m and so on

        def _par_or_zeros(par):
            # Return the parameter if it is not None, otherwise zeros
            return par if par is not None else (0.0, 0.0, 0.0)

        self.mass_plus = mass_plus
        self.mass_minus = mass_minus
        self.coordinate_distance = coordinate_distance

        self.par_b = coordinate_distance / 2

        total_mass = self.mass_plus + self.mass_minus

        self.coord_x_plus = (
            self.coordinate_distance * self.mass_minus / total_mass
        )
        self.coord_x_minus = (
            -self.coordinate_distance * self.mass_plus / total_mass
        )

        self.momenta_plus = _par_or_zeros(momenta_plus)
        self.momenta_minus = _par_or_zeros(momenta_minus)
        self.chi_plus = _par_or_zeros(chi_plus)
        self.chi_minus = _par_or_zeros(chi_minus)

        self.S_plus = tuple(cp * self.mass_plus ** 2 for cp in self.chi_plus)
        self.S_minus = tuple(
            cm * self.mass_minus ** 2 for cm in self.chi_minus
        )

        if center_offset is None:
            self.center_offset = (self.coord_x_plus - self.par_b, 0.0, 0.0)
        else:
            self.center_offset = center_offset

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        def assign_parameter(
            param, value, which_bh=None, vector=False, thorn=None
        ):
            if thorn is None:
                thorn = "TwoPunctures"

            if which_bh is None:
                return f"{thorn}::{param} = {value}"

            # which_bh is either _plus or _minus
            if vector is None:
                return f"{thorn}::{param}{which_bh} = {value}"

            if vector:
                return "\n".join(
                    [
                        f"{thorn}::{param}{which_bh}[{index}] = {value[index]}"
                        for index in range(3)
                    ]
                )
            return f"TwoPunctures::{param}{which_bh} = {value}"

        ret = []

        ret.append(
            assign_parameter("initial_data", "twopunctures", thorn="ADMBase")
        )
        ret.append(
            assign_parameter(
                "initial_lapse", "twopunctures-averaged", thorn="ADMBase"
            )
        )
        ret.append(assign_parameter("initial_shift", "zero", thorn="ADMBase"))
        ret.append(
            assign_parameter("initial_dtlapse", "zero", thorn="ADMBase")
        )
        ret.append(
            assign_parameter("initial_dtshift", "zero", thorn="ADMBase")
        )

        ret.append(assign_parameter("give_bare_mass", "no"))
        ret.append(assign_parameter("par_b", self.par_b))
        ret.append(
            assign_parameter("target_m", self.mass_plus, "_plus", False)
        )
        ret.append(
            assign_parameter("target_m", self.mass_minus, "_minus", False)
        )
        ret.append(assign_parameter("par_P", self.momenta_plus, "_plus", True))
        ret.append(
            assign_parameter("par_P", self.momenta_minus, "_minus", True)
        )
        ret.append(assign_parameter("par_S", self.S_plus, "_plus", True))
        ret.append(assign_parameter("par_S", self.S_minus, "_minus", True))
        ret.append(
            assign_parameter("center_offset", self.center_offset, "", True)
        )

        return "\n".join(ret)
