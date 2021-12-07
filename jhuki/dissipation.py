#!/usr/bin/env python3

# Copyright (C) 2020-2021 Gabriele Bozzola
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


"""The :py:mod:`~.dissipation` module takes care of setting parameters for the
dissipation thorn. The main class is :py:class:`~.Dissipation`, but it is more
convenient to use helper functions to setup this object. For example, to set
up dissipation for a simulation with Lean, the simplest way is:

.. code-block:: python

   dis = add_gauge(add_Lean(create_dissipation_from_grid(grid, 0.3)))

where 0.3 is the value of ``eps_dis`` at the finest level.

"""

import enum
from functools import lru_cache

from jhuki.base import BaseThorn


class Dissipation(BaseThorn):
    """Class that describes the dissipation. Best used in combination with helper
    functions.

    This class is immutable.

    :ivar epsdis_per_level: Dictionary that maps the value of the dissipation on
                            each refinement level.
    :type epsdis_per_level: dict
    :ivar order: Order of the dissipation.
    :type order: int
    :ivar variables: List of variables to which the dissipation has to be applied.
    :type variables: list of str

    """

    def __init__(self, epsdis_per_level, order, variables):
        """It is your responsibility to provide meaningful eps_dis_per_level.

        :param epsdis_per_level: Dictionary that maps the value of the
                                 dissipation on each refinement level.
        :type epsdis_per_level: dict
        :param order: Order of the dissipation.
        :type order: int
        :param variables: List of variables to which the dissipation has to be
                          applied.
        :type variables: list of str

        """
        self.epsdis_per_level = epsdis_per_level
        self.order = order
        self.variables = variables

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        eps_dis_str = "\n".join(
            f"Dissipation::epsdis_for_level[{level}] = {epsdis}"
            for level, epsdis in self.epsdis_per_level.items()
        )

        vars_str = "\n".join(self.variables)

        return f"""\
Dissipation::order = {self.order}
Dissipation::vars  = "{vars_str}"
{eps_dis_str}"""


class DissPrescription(enum.Enum):
    """DissPrescription for how to set up the dissipation.

    Three different prescriptions are available:
    - ``const``: set the same ``eps_dis`` everywhere,
    - ``dtfac``: set ``eps_dis`` proportional to the local dtfac,
    - ``continuous``: set ``eps_dis`` proportional to the local grid spacing
                      to the power of the order

    """

    const = enum.auto()
    dtfact = enum.auto()
    continuous = enum.auto()


def create_dissipation_from_grid(
    grid, eps_dis_finest, variables=None, prescription=DissPrescription.dtfact
):
    """Create a :py:class:`~.Dissipation` object from a given :py:class:`~.Grid`.


    The order of the dissipation is deduced from the number of ghost zones.

    :param grid: Simulation grid.
    :type grid: :py:class:`~.Grid`
    :param epsdis_per_level: Dictionary that maps the value of the dissipation
                             on each refinement level.
    :type epsdis_per_level: dict
    :param order: Order of the dissipation.
    :type order: int
    :param variables: List of variables to which the dissipation has to be
                      applied.
    :type variables: list of str
    :param prescription: How to set ``eps_dis``.
    :type prescription: :py:class:`~.DissPrescription`

    :returns: Dissipation object.
    :rtype: :py:class:`~.Dissipation`

    """
    if not isinstance(prescription, DissPrescription):
        raise ValueError(
            f"Unknown prescription {prescription}. "
            "Accepted values: {list(DissPrescription.__members__.keys())}"
        )

    order = 2 * grid.num_ghost - 1
    variables = [] if variables is None else variables

    # If we have N levels, we need to set up dissipation on N + 1 zones
    num_levels = grid.max_num_refinement_levels + 1

    if prescription == DissPrescription.const:
        # Same everywhere
        eps_dis_levels = {level: eps_dis_finest for level in range(num_levels)}
    elif prescription == DissPrescription.dtfact:
        # Same everywhere, halved by two every time a level is synced
        # For example, if have 4 levels, 2 synced, then we should have
        # [1/4, 1/2, 1, 1] * eps_dis

        num_levels_synced = grid.num_levels_with_dt_coarse

        # First, the ones with the values halved
        eps_dis_levels = {
            level: eps_dis_finest / 2 ** (num_levels_synced - level)
            for level in range(num_levels_synced)
        }

        # Next, we add all the ones that are eps_dis_finest
        eps_dis_levels.update(
            {
                level: eps_dis_finest
                for level in range(
                    num_levels_synced,
                    num_levels,
                )
            }
        )

    elif prescription == DissPrescription.continuous:
        # Proportional to the local size, so there's a factor of 2**order for every
        # refinement level
        eps_dis_levels = {
            level: eps_dis_finest / (2 ** (num_levels - level - 1)) ** order
            for level in range(num_levels)
        }

    return Dissipation(eps_dis_levels, order, variables)


def add_to_dissipation(variables, doc):
    """Return a function that adds variables to the given Dissipation."""

    def inner(dissipation):
        __doc__ = doc  # noqa: F841
        return Dissipation(
            dissipation.epsdis_per_level,
            dissipation.order,
            dissipation.variables + variables,
        )

    return inner


_lean_vars = [
    "LeanBSSNMoL::conf_fac",
    "LeanBSSNMoL::hmetric",
    "LeanBSSNMoL::hcurv",
    "LeanBSSNMoL::trk",
    "LeanBSSNMoL::gammat",
]

_gauge_vars = ["ADMBase::lapse", "ADMBase::shift"]

_proca_vars = [
    "ProcaBase::Ei",
    "ProcaBase::Ai",
    "ProcaBase::Aphi",
    "ProcaBase::Zeta",
]

add_Lean = add_to_dissipation(
    _lean_vars,
    """Return a new :py:class:`~.Dissipation` with Lean's variables added.

    :param dissipation: Dissipation that has to be used as base.
    :type dissipation: :py:class:`~.Dissipation`
    :returns: Dissipation with Lean's variables added.
    :rtype: :py:class:`~.Dissipation`""",
)

add_gauge = add_to_dissipation(
    _gauge_vars,
    """Return a new :py:class:`~.Dissipation` with guage variables added.

    :param dissipation: Dissipation that has to be used as base.
    :type dissipation: :py:class:`~.Dissipation`
    :returns: Dissipation with guage variables added.
    :rtype: :py:class:`~.Dissipation`""",
)

add_Proca = add_to_dissipation(
    _proca_vars,
    """Return a new :py:class:`~.Dissipation` with Proca variables added.

    :param dissipation: Dissipation that has to be used as base.
    :type dissipation: :py:class:`~.Dissipation`
    :returns: Dissipation with Proca variables added.
    :rtype: :py:class:`~.Dissipation`""",
)
