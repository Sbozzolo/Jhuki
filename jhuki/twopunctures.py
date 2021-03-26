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


def prepare_quasicircular_inspiral(
    mass_ratio,
    coordinate_distance,
    total_bare_mass=1,
    chi_plus=(0, 0, 0),
    chi_minus=(0, 0, 0),
):
    """Return a :py:class:`~.TwoPunctures` that describes a quasi-circular inspiral.
    We always assume that the plus puncture is the most massive one.

    :param mass_ratio: Mass ratio of the binary.
    :type mass_ratio: float
    :param coordinate_distance: Initial coordinate separation.
    :type coordinate_distance: float
    :param total_bare_mass: Rescale masses and momenta so that the total bare mass of the
                            system is this.
    :type total_bare_mass: float

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

    Pt, Pr = Pt * total_bare_mass, Pr * total_bare_mass

    momenta_plus = (-Pr, Pt, 0)
    momenta_minus = (Pr, -Pt, 0)

    mass_plus = mass_ratio / (1 + mass_ratio) * total_bare_mass
    mass_minus = 1 / (1 + mass_ratio) * total_bare_mass

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

    :ivar swap_xz: If True, activate the ``swap_xz`` parameter in TwoPunctures.
    :vartype swap_xz: bool
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
        swap_xz=False,
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

        :param swap_xz: If True, activate the ``swap_xz`` parameter in TwoPunctures.
        :type swap_xz: bool
        """

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

        self.center_offset = (self.coord_x_plus - self.par_b, 0.0, 0.0)

        self.swap_xz = swap_xz

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        def assign_parameter(param, value, which_bh=None):
            if which_bh is None:
                return f"TwoPunctures::{param} = {value}"

            # which_bh is either _plus or _minus
            return "\n".join(
                [
                    f"TwoPunctures::{param}{which_bh}[{index}] = {value[index]}"
                    for index in range(3)
                ]
            )

        ret = []

        ret.append(
            f"""\
ADMBase::initial_data = "twopunctures"
ADMBase::initial_lapse = "twopunctures-averaged"
ADMBase::initial_shift = "zero"
ADMBase::initial_dtlapse = "zero"
ADMBase::initial_dtshift = "zero"

TwoPunctures::give_bare_mass = "no"
TwoPunctures::par_b = {self.par_b}
TwoPunctures::target_m_plus = {self.mass_plus}
TwoPunctures::target_m_minus = {self.mass_minus}
TwoPunctures::par_m_plus = {self.mass_plus}
TwoPunctures::par_m_minus = {self.mass_minus}"""
        )

        ret.append(assign_parameter("par_P", self.momenta_plus, "_plus"))
        ret.append(assign_parameter("par_P", self.momenta_minus, "_minus"))
        ret.append(assign_parameter("par_S", self.S_plus, "_plus"))
        ret.append(assign_parameter("par_S", self.S_minus, "_minus"))
        ret.append(assign_parameter("center_offset", self.center_offset, ""))

        if self.swap_xz:
            ret.append(assign_parameter("swap_xz", "yes"))

        return "\n".join(ret)
