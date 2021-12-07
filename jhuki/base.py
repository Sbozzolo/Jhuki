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

"""The :py:mod:`~.base` module contains basic functions and classes used
everywhere in the codebase.


"""

from abc import ABC, abstractmethod


class BaseThorn(ABC):
    """Base abstract class for thorn definitions.

    It defines the ``__str__`` method from the ``parfile_code`` one.

    Derived classes have to define a ``parfile_code`` method

    """

    @abstractmethod
    def parfile_code(self):  # pragma: no cover
        raise NotImplementedError

    def __str__(self):
        return self.parfile_code
