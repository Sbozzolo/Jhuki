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

from jhuki.twopunctures import TwoPunctures


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
        center_offset=None,
        charge_plus=None,
        charge_minus=None,
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

        super().__init__(
            mass_plus,
            mass_minus,
            coordinate_distance,
            momenta_plus=momenta_plus,
            momenta_minus=momenta_minus,
            chi_plus=chi_plus,
            chi_minus=chi_minus,
            center_offset=center_offset,
        )

        self.charge_plus = charge_plus if charge_plus is not None else 0
        self.charge_minus = charge_minus if charge_minus is not None else 0

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        uncharged_parfile = super().parfile_code.replace(
            "TwoPunctures", "TwoChargedPunctures"
        )

        charged_parfile = (
            f"{uncharged_parfile}\n"
            f"TwoChargedPunctures::par_q_plus = {self.charge_plus}\n"
            f"TwoChargedPunctures::par__plus = {self.charge_plus}"
        )

        return charged_parfile
