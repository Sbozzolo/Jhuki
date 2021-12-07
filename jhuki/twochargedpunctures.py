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

"""The :py:mod:`~.twochargedpunctures` module extends :py:mod:`~.twopunctures`
 for the case of charged black holes. The interface is nearly identical.
"""

from functools import lru_cache
from math import sqrt

from jhuki.externals.nrpypn import compute_quasicircular_momenta
from jhuki.twopunctures import TwoPunctures


def prepare_quasicircular_inspiral(
    mass_ratio,
    coordinate_distance,
    total_bare_mass=1,
    chi_plus=(0, 0, 0),
    chi_minus=(0, 0, 0),
    lambda_plus=0,
    lambda_minus=0,
    **kwargs,
):
    """Return a :py:class:`~.TwoChargedPunctures` that describes a quasi-circular inspiral.
    We always assume that the plus puncture is the most massive one.

    We compute the momenta starting from the uncharged case and scaling by
    sqrt(1 - lambda1 * lambda2). This is essentially 0 PN. It works well for low
    charge, but not that much for higher charge.

    :param mass_ratio: Mass ratio of the binary.
    :type mass_ratio: float
    :param lambda_plus: Charge-to-mass-ratio of the black hole on the positive
                        side of the x (or z) axis.
    :type lambda_plus: float
    :param lambda_minus: Charge-to-mass-ratio of the black hole on the negative
                         side of the x (or z) axis.
    :type lambda_minus: float
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

    Unknown arguments are passed to :py:class:`~.TwoPunctures`.

    :returns: A :py:class:`~.TwoPunctures` for a quasi-circular inspiral.
    :rtype: :py:class:`~.TwoPunctures`

    """

    if mass_ratio < 1:
        mass_ratio = 1 / mass_ratio

    Pt, Pr = compute_quasicircular_momenta(
        mass_ratio, coordinate_distance, chi_plus, chi_minus
    )

    mass_plus = mass_ratio / (1 + mass_ratio) * total_bare_mass
    mass_minus = 1 / (1 + mass_ratio) * total_bare_mass

    factor = sqrt(1 - lambda_plus * lambda_minus)

    Pt, Pr = Pt * total_bare_mass * factor, Pr * total_bare_mass * factor

    momenta_plus = (-Pr, Pt, 0)
    momenta_minus = (Pr, -Pt, 0)

    return TwoChargedPunctures(
        mass_plus,
        mass_minus,
        coordinate_distance,
        momenta_plus=momenta_plus,
        momenta_minus=momenta_minus,
        chi_plus=chi_plus,
        chi_minus=chi_minus,
        charge_plus=lambda_plus * mass_plus,
        charge_minus=lambda_minus * mass_minus,
        give_bare_mass=False,
        **kwargs,
    )


class TwoChargedPunctures(TwoPunctures):
    """The :py:class:`~.TwoChargedPunctures` extends :py:class:`~.TwoPunctures` for
    the case charged black holes (to use with the ``TwoChargedPunctures`` thorn).

    As :py:class:`~.TwoPunctures`, this class is immutable.

    :ivar charge_plus: Charge of the black hole on the positive side of the x (or z) axis.
    :vartype charge_plus: float
    :ivar charge_minus: Charge of the black hole on the negative side of the x (or z) axis.
    :vartype charge_minus: float
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
        charge_plus=0,
        charge_minus=0,
        swap_xz=False,
        give_bare_mass=False,
    ):
        """Constructor.

        :param mass_plus: Mass of the black hole on the positive side of the x (or z) axis.
        :type mass_plus: float
        :param mass_minus: Mass of the black hole on the negative side of the x (or z) axis.
        :type mass_minus: float

        :param charge_plus: Charge of the black hole on the positive side of the x (or z) axis.
        :type charge_plus: float
        :param charge_minus: Charge of the black hole on the negative side of the x (or z) axis.
        :type charge_minus: float

        :param coordinate_separation: Initial coordinate separation.
        :type coordinate_separation: float
        :param par_b: Half of the initial separation.
        :type par_b: float

        :param momenta_plus: Array with the linear momenta along the three directions for the
                            black hole on the positive side of the x (or z) axis.
                             If ``swap_xz`` is True, these are the values after swapping.
                             For example ``(0, 0, 0.5)`` is a vector along the z axis.
                             (This is the opposite of what happens in ``TwoPunctures``, where
                              you have to define the values before swapping)
        :type momenta_plus: list/tuple
        :param momenta_minus: Array with the linear momenta along the three directions for the
                             black hole on the positive side of the x (or z) axis.
                             If ``swap_xz`` is True, these are the values after swapping.
                             For example ``(0, 0, 0.5)`` is a vector along the z axis.
                             (This is the opposite of what happens in ``TwoPunctures``, where
                              you have to define the values before swapping)
        :type momenta_minus: list/tuple

        :param chi_plus: Array with the dimensionless_spin along the three directions for the
                        black hole on the positive side of the x (or z) axis.
                             If ``swap_xz`` is True, these are the values after swapping.
                             For example ``(0, 0, 0.5)`` is a vector along the z axis.
                             (This is the opposite of what happens in ``TwoPunctures``, where
                              you have to define the values before swapping)
        :type chi_plus: list/tuple
        :param chi_minus: Array with the dimensionless spin along the three directions for the
                         black hole on the positive side of the x (or z) axis.
                             If ``swap_xz`` is True, these are the values after swapping.
                             For example ``(0, 0, 0.5)`` is a vector along the z axis.
                             (This is the opposite of what happens in ``TwoPunctures``, where
                              you have to define the values before swapping)
        :type chi_minus: list/tuple

        :param swap_xz: If True, activate the ``swap_xz`` parameter in TwoPunctures.
        :type swap_xz: bool

        :param give_bare_mass: If True, set this parameter to True in the parfile.
        :type give_bare_mass: bool
        """

        super().__init__(
            mass_plus,
            mass_minus,
            coordinate_distance,
            momenta_plus=momenta_plus,
            momenta_minus=momenta_minus,
            chi_plus=chi_plus,
            chi_minus=chi_minus,
            swap_xz=swap_xz,
            give_bare_mass=give_bare_mass,
            initial_alpha="psi^n",
        )

        self.charge_plus = charge_plus
        self.charge_minus = charge_minus

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        uncharged_parfile = (
            super()
            .parfile_code.replace("TwoPunctures", "TwoChargedPunctures")
            .replace("twopunctures", "twochargedpunctures")
        )

        charged_parfile = (
            f"{uncharged_parfile}\n"
            f"TwoChargedPunctures::par_q_plus = {self.charge_plus}\n"
            f"TwoChargedPunctures::par_q_minus = {self.charge_minus}"
        )

        return charged_parfile
