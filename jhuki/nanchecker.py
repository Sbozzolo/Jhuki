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


"""The :py:mod:`~.nanchecker` module sets the parameter for NaNChecker.

This is a simple module and the only non-trivial function is
:py:func:`~.create_nanchecker_from_grid`, which automatically sets the
``check_every`` parameter.

"""

from functools import lru_cache

from jhuki.base import BaseThorn


class NaNChecker(BaseThorn):
    """Class that describes NaNChecker.

    This class is immutable.

    :ivar action_if_found: What to do if a NaN is found.
    :type action_if_found: str
    :ivar check_every: Check for NaNs every N iterations.
    :type check_every: int
    :ivar variables: List of variables for which NaNs have to be checked.
    :type variables: list of str

    """

    def __init__(self, variables, check_every, action_if_found="terminate"):
        """Constructor.

        :param check_every: Check for NaNs every N iterations.
        :type check_every: int
        :param variables: List of variables for which NaNs have to be checked.
        :type variables: list of str
        :param action_if_found: What to do if a NaN is found.
        :type action_if_found: str

        """
        self.variables = variables
        self.check_every = check_every
        self.action_if_found = action_if_found

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        vars_str = "\n".join(self.variables)

        return f"""\
NaNChecker::check_every = {self.check_every}
NaNChecker::action_if_found = {self.action_if_found}
NaNChecker::check_vars  = "{vars_str}"\
"""


def create_nanchecker_from_grid(grid, variables, action_if_found="terminate"):
    """Create a :py:class:`~.NaNChecker` object from a given :py:class:`~.Grid`.

    The only information used is when refinement levels are synced.

    :param grid: Simulation grid.
    :type grid: :py:class:`~.Grid`
    :param variables: List of variables for which NaNs have to be checked.
    :type variables: list of str
    :param action_if_found: What to do if a NaN is found.
    :type action_if_found: str

    :returns: NaNChecker object.
    :rtype: :py:class:`~.NaNChecker`

    """
    return NaNChecker(
        variables, grid.rl_synced_every, action_if_found=action_if_found
    )
