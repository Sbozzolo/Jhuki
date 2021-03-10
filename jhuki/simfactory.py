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

"""The :py:mod:`~.simfactory` module provides functions to interface with
SimFactory.

The available functions are: - :py:func:`~.simfactory_option`, which helps
setting up options that are controllable by SimFactory.

"""


def simfactory_option(name, default):
    """Set a variable to a SimFactory-controllable parameter.

    SimFactory can define arbitrary parameters using the ``--define`` keyword.
    With this function, you can define a variable that takes a specified default
    value, unless the keyword is used.

    For example:

    .. code-block::python

        par_b = simfactory_option("@PARB@", 6.5)

    par_b is a Python variable with the value 6.5 unless ``--define PARB=XX`` is
    passed to SimFactory. The way this works is that SimFactory performs a
    substitution of objects with the double @ symbols. For this reason, the
    name parameter to this function has to be surrounded by @.

    :param name: Name of the SimFactory define key. It has to be with surrounded
                 by @, for example, @XMAX@.
    :type name: str
    :param default: Default value that the variable has to have if define is not
                    specified.
    :type default: int, float, or str

    :returns: The default value if the define key is not passed to SimFactory,
              otherwise the value defined with SimFactory.
    :rtype: float or str (integers are converted )
    """
    if not isinstance(name, str):
        raise TypeError(f"name has to a string, it is {name}")

    # If default is a string, we will keep it a string. If it is an int, we
    # convert it to float. Often one sets parameter like par_b = 6, but here 6
    # has to be intended as a float, or it can cause troubles for other
    # calculations.
    type_ = float if isinstance(default, int) else type(default)

    return type_(name) if name[0] != "@" else type_(default)
