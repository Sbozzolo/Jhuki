#!/usr/bin/env python3

# Copyright (C) 2020-2022 Gabriele Bozzola
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


"""The :py:mod:`~.grid` module provides infrastructure to work with the grid
structure of a simulation.

The module provides a class :py:class:`~.RefinementCenter` which contains all
information for the refinement structure around a given refinement center. One
or more :py:class:`~.RefinementCenter` objects form a :py:class:`~.Grid`, which
has the entire information on a computational grid. The object can produce
parfile code or can be printed to give an overview of the grid structure.

One of the most important functions provided in this module is
:py:func:`~.set_dt_max_grid`. The function takes two arguments: a
:py:class:`~.Grid` and the maximum timestep that you want the grid to have. The
return value is a new :py:class:`~.Grid` which is identical to the input one,
with the exception that the timestep of the outermost refinement levels is
synced so that the maximum timestep on the grid is the one you provided. This
function is important because it can be used to produce grids that satisfy CFL
conditions (like the ones due to damp_lorenz or eta_beta).

Another useful function is :py:func:`~.create_twopunctures_grid`, which takes a
:py:class:`~.TwoPunctures` and returns a grid built on top of it.

Warning: at the moment we only support grids in which the three dimensions are
the same in resolution and extent.

"""

import logging
from functools import lru_cache
from math import ceil, sqrt
from typing import Iterable, List, Optional, Union

from jhuki.base import BaseThorn
from jhuki.twochargedpunctures import TwoChargedPunctures
from jhuki.twopunctures import TwoPunctures


def create_twopunctures_grid(
    two_punctures: Union[TwoPunctures, TwoChargedPunctures],
    points_on_horizon_radius: int,
    minimum_outer_boundary: float,
    cfl_num: float = 0.45,
    horizon_padding_points: int = 5,
    skip_radii: Optional[Iterable[int]] = None,
    **kwargs,
):
    """Create a grid centered around a given two punctures.

    The grid has three identical refinement centers: two on the punctures and
    one in the center.

    The construction of the grid starts by estimating the horizon radius. Then,
    the horizon is covered with ``points_on_horizon_radius``. From there,
    refinement levels are added until the ``minimum_outer_boundary`` is reached.
    Radii are placed at powers of 2 from the central one, unless the index is
    in the ``skip_radii`` list. For example, if ``skip_radii = [2, 3]``, and the
    refinement boundaries are supposed to be ``[1, 2, 4, 8, 16]``, the actual
    one will be ``[1, 2, 16]``. This is useful to encompass regions of space with
    the same resolution. Indices are counted from the innermost.

    This function does nothing fancy about the mass ratio: it simply considers
    the smaller horizon and creates the grid starting from there.

    Unknown arguments are passed to the constructor of :py:class:`~.Grid`.

    :param two_punctures: Two punctures for which the grid has to be created.
    :type two_punctures: :py:class:`~.TwoPunctures` or
                         :py:class:`~.TwoChargedPunctures`
    :param points_on_horizon_radius: Desired number of grid points on the
                                     radius of the smaller horizon.
    :type points_on_horizon_radius: int
    :param minimum_outer_boundary: Minimum outer boundary that the grid has
                                   to have
    :type minimum_outer_boundary: float
    :param cfl_num: CFL factor for the finest refinement level
    :type cfl_num: float

    :param horizon_padding_points: Add this number of points to the innermost
                                   refinement level.
    :type horizon_padding_points: int

    :param skip_radii: Do not add to the grid refinement boundaries with these
                       indices.
    :type skip_radii: list of int

    :returns: Grid built on top of the given two punctures.
    :rtype: :py:class:`~.Grid`

    """
    # First, we estimate the radius of the horizons

    mp = two_punctures.mass_plus
    mm = two_punctures.mass_minus

    # We take 1e-10 as a small number below which the punctures are assumed to
    # be massless
    mp_is_nonzero = mp > 1e-10
    mm_is_nonzero = mm > 1e-10

    if isinstance(two_punctures, TwoChargedPunctures):
        qp = two_punctures.charge_plus
        qm = two_punctures.charge_minus
    else:
        qp, qm = 0, 0

    # Spin
    ap_sq = sum(a_i**2 for a_i in two_punctures.a_plus)
    am_sq = sum(a_i**2 for a_i in two_punctures.a_minus)

    # For the size at t = 0, we would divide this by 2. In practice, the
    # horizons expand in a few Ms to reach the guess_r size, so that's the more
    # relevant quantity.

    # This is used only to compute dx, so we set the radius to a huge number
    # when the mass is zero. The `min` below will select the correct puncture.
    guess_r_plus = sqrt(mp**2 - ap_sq - qp**2) if mp_is_nonzero else 1e100
    guess_r_minus = sqrt(mm**2 - am_sq - qm**2) if mm_is_nonzero else 1e100

    # Now, we can find the resolution needed to cover the smaller horizon with
    # the desired number of points

    # We truncate the dx to four decimals so that it is easy to ensure that we
    # have an integer number of points
    dx_fine = round(
        min(guess_r_plus, guess_r_minus) / points_on_horizon_radius, 4
    )

    # radius0 is the radius of the innermost refinement level, the one that
    # encompasses the horizon. We add 5 points so that we can ensure that the
    # horizon is well contained
    radius0 = dx_fine * (points_on_horizon_radius + horizon_padding_points)

    # We can build up the hierarchy of levels. We add one level at the time
    # until we are at the desired outer boundary.
    radius_number = 0
    radius = radius0
    # Here is were the save the various radii
    radii = []

    if skip_radii is None:
        skip_radii = []

    while radius < minimum_outer_boundary:
        radius = radius0 * 2**radius_number
        if radius_number not in skip_radii:
            radii.append(radius)
        radius_number += 1

    # We do not want to include the outer boundary
    radii, outer_boundary = radii[:-1], radii[-1]

    if two_punctures.swap_xz:
        pp = (0, 0, two_punctures.coord_x_plus)
        pm = (0, 0, two_punctures.coord_x_minus)
    else:
        pp = (two_punctures.coord_x_plus, 0, 0)
        pm = (two_punctures.coord_x_minus, 0, 0)

    centers: List[RefinementCenter] = []

    if mp_is_nonzero:
        centers.append(
            RefinementCenter(
                radii,
                dx_fine,
                cfl_num,
                center_num=len(centers) + 1,
                position=pp,
            )
        )
    if mm_is_nonzero:
        centers.append(
            RefinementCenter(
                radii,
                dx_fine,
                cfl_num,
                center_num=len(centers) + 1,
                position=pm,
            )
        )

    if mp_is_nonzero and mm_is_nonzero:

        # The refinement level in the center is needed mostly to ensure that the
        # interaction is well resolved. We define "well resolved" as in the coarsest
        # resolution that we allow is dx = dist / 80, meaning, we have at least 80
        # points to cover the distance.
        min_dx_on_separatrix = two_punctures.coordinate_distance / 80

        # Now, we find what is the dx and number of the first level that satisfied
        # this condition

        target_index = None

        # Note, center1.dx is boundary-to-center, which is the opposite direction of
        # everything else, so, we reverse it. Here, we traverse from center to
        # boundary and we find the first level for which the resolution is not good
        # enough.
        dxs_center_to_boundary = list(reversed(centers[0].dx))

        for ref_index, dx in enumerate(dxs_center_to_boundary):
            if dx > min_dx_on_separatrix:
                # We subtract one because the last good level was the previous
                target_index = ref_index - 1
                break

        if target_index is None:
            raise RuntimeError(
                "There is not enough points to resolve the separatrix well"
            )

        dx_fine_center = dxs_center_to_boundary[target_index]

        centers.append(
            RefinementCenter(
                radii[target_index:],
                dx_fine_center,
                cfl_num,
                center_num=len(centers) + 1,
                position=(0, 0, 0),
            )
        )

    # We adjust the outer boundary to ensure that it is multiple of dx_coarse
    # This gives the first larger value that is multiple of dx_coarse
    dx_coarse = centers[0].dx_coarse
    outer_boundary = ceil(outer_boundary / dx_coarse) * dx_coarse

    return Grid(centers, outer_boundary, **kwargs)


def set_dt_max_grid(grid, dt_max):
    """Return a new :py:class:`~.Grid` with maximum timestep smaller than
    the given one.

    This is achieved by applying :py:func:`~.set_dt_max` to each refinement
    center.

    This function is very important! Several equations have a maximum timestep
    allowed for numerical stability (e.g., the shift gauge evolution equation
    and its parameter eta_beta).

    :param grid: Grid that has to be adjusted.
    :type grid: :py:class:`~.Grid`
    :param dt_max: Maximum timestep allowed.
    :type dt_max: float

    :returns: Grid with maximum timestep smaller than the given one.
    :rtype: :py:class:`~.Grid`
    """

    # TODO: Add type checks

    # When we work with grids where different refinement centers have a
    # different number of levels, it can happen that some centers end up
    # completely synced up. But at the end of the day, we only care about the
    # center with the largest number of levels. we have already checked that all
    # the refinement centers are compatible. So here we can just focus on the
    # one with the most number of levels.

    new_center_with_most_levels = set_dt_max(
        grid.center_with_most_levels, dt_max
    )
    num_levels_synced = new_center_with_most_levels.num_levels_with_dt_coarse

    new_refinement_centers = set()

    # TODO: FINISH HERE

    for ref_center in grid.refinement_centers:
        # We set dx and cfl reading it from new_center_with_most_levels.
        # The -1 is because we start counting from zero
        index_level = ref_center.num_refinement_levels - 1

        # Here we take care of the two options: we have to sync all the levels,
        # or we have to sync only some of the levels
        num_levels_with_dt_coarse = (
            index_level + 1
            if index_level < num_levels_synced
            else num_levels_synced
        )

        new_ref_center = RefinementCenter(
            refinement_radii=ref_center.refinement_radii,
            dx_fine=new_center_with_most_levels.dx[index_level],
            cfl_fine=new_center_with_most_levels.cfl[index_level],
            center_num=ref_center.center_num,
            position=ref_center.position,
            num_levels_with_dt_coarse=num_levels_with_dt_coarse,
        )

        new_refinement_centers.add(new_ref_center)

    new_refinement_centers = tuple(new_refinement_centers)
    return Grid(
        new_refinement_centers,
        outer_boundary=grid.outer_boundary,
        num_ghost=grid.num_ghost,
        symmetry=grid.symmetry,
        tiny_shift=grid.tiny_shift,
    )


def set_dt_max(ref_center, dt_max):
    """Return a new :py:class:`~.RefinementCenter` with maximum timestep smaller than
    the given one. This is achieved by changing the parameter
    `~.num_levels_with_dt_coarse`.

    This function is smart enough to preserve existing synced levels, so you can
    apply it as many times as you want (for diverse conditions).

    This function is very important! Several equations have a maximum timestep
    allowed for numerical stability (e.g., the shift gauge evolution equation
    and its parameter eta_beta).

    :param ref_center: Grid that has to be adjusted.
    :type ref_center: :py:class:`~.RefinementCenter`
    :param dt_max: Maximum allowed timestep.
    :type dt_max: float

    :returns: New grid with maximum timestep smaller than ``dt_max``.
    :rtype: :py:class:`~.RefinementCenter`

    """
    # TODO: Add type checks

    logging.debug("Working with grid:")
    logging.debug(ref_center)

    # What level has dt < dt_max?
    #
    # ref_center.dt is from the coarsest, so if we can walk through the list and
    # see when dt_max >= dt. In this, we have to pay attention to the fact that
    # ref_center might have num_levels_with_dt_coarse != 1. In this case, we have
    # to skip the refinement levels that are synced and starting from the latest
    # synced one. This explains the `ref_center.num_levels_with_dt_coarse - 1`, if
    # we have three levels evolving with dt_coarse, the index of the last level
    # of this group is 2.
    for level in range(
        ref_center.num_levels_with_dt_coarse - 1,
        ref_center.num_refinement_levels,
    ):
        dt = ref_center.dt[level]
        if dt_max >= dt:
            # Let us go through an example. Say we have 3 refinement levels (2
            # refinement radii), with dx = [4, 2, 1] and dt = [2, 1, 0.5]. We
            # then set dt_max = 0.75, what we want is that all the three levels
            # evolve with the same dt = 0.5. num_levels_with_dt_coarse has to
            # be 3. This explains the `+ 1`; we always have to at least one
            # level evolving at dt_coarse
            num_levels_with_dt_coarse = level + 1
            # Here we create a new object and what we change is only
            # num_levels_with_dt_coarse
            logging.debug(
                f"New num_levels_with_dt_coarse: {num_levels_with_dt_coarse}"
            )
            return RefinementCenter(
                refinement_radii=ref_center.refinement_radii,
                dx_fine=ref_center.dx_fine,
                cfl_fine=ref_center.cfl_fine,
                center_num=ref_center.center_num,
                position=ref_center.position,
                num_levels_with_dt_coarse=num_levels_with_dt_coarse,
            )

    raise ValueError(
        f"Given dt_max {dt_max} is larger than the smallest timestep in "
        f"the grid {ref_center.dt_fine}"
    )


class RefinementCenter(BaseThorn):
    """The :py:class:`~.RefinementCenter` class contains all the information related
    to the grid structure around one single refinement center. This is not the
    entire grid structure for a simulation (unless you have only one refinement
    center).

    This class is immutable.

    To understand this class always keep in mind the following example, where we
    draw a simple grid with three refinement levels (two refinement radii). The
    numbers 500, 250, 125, 0 are hypothetical refinement radii that you would
    put in CarpetRegrid.

    ..code-block::

        |------------------|------------------|------------------|
       500                250                125                 0
       O.B.             R.R. 1             R.R. 2
              R.L. 0              R.L. 1             R.L. 2


    :ivar refinement_radii: Radii over which a change in refinement happens.
                            This is what you would write in CarpetRegrid.
                            Refinement levels are always stored from the largest
                            to the smallest (as in Carpet).
    :vartype refinement_radii: list of float

    :ivar center_num: Identifier of the current mesh refinement hierarchy in the
                      context of the entire grid. In CarpetRegrid, you specify
                      num_centers, then each hierarchy is identified with
                      numbers from 1 to num_centers. This is that number.
    :vartype center_num: int

    :ivar position: (x, y, z) position of the refinement center in the overall
                     grid.
    :vartype position: list or tuple of float

    :ivar dx_fine: Resolution at the finest level (near the center).
                   We assume that the three directions have the same resolution.
    :vartype dx_fine: float

    :ivar cfl_fine: CFL number at the finest level (near the center).
                      This is dt_finest/dx_finest.
    :vartype cfl_fine: float

    :ivar dt_fine: Timestep at the finest level (near the center).
    :vartype dt_fine: float

    :ivar num_levels_with_dt_coarse: How many refinement levels should have the
                                     same dt as the coarsest refinement level.
                                     This will translate to
                                     Carpet::time_refinement_factors. The
                                     minimum is 1 (the coarsest level itself).
    :vartype num_levels_with_dt_coarse: int

    :ivar num_refinement_radii: Number of radii over which a change in
                                refinement happens. This is not the number of
                                refinement levels (while the two are off by
                                one). This is the number of radii in
                                CarpetRegrid.
    :vartype num_refinement_radii: int

    :ivar num_refinement_levels: Number of distinct refinement levels. This is
                                 the same as the number of different resolutions
                                 over a grid and it is equal to the number of
                                 refinement radii + 1.

    :vartype num_refinement_levels: int

    """

    def __init__(
        self,
        refinement_radii,
        dx_fine,
        cfl_fine,
        center_num=1,
        position=(0, 0, 0),
        num_levels_with_dt_coarse=1,
    ):
        """Constructor.

        :param refinement_radii: Radii over which a change in refinement happens.
                                 This is what you would write in CarpetRegrid.
        :type refinement_radii: list of float

        :param center_num: Identifier of the current mesh refinement hierarchy
                           in the context of the entire grid. In CarpetRegrid,
                           you specify num_centers, then each hierarchy is
                           identified with numbers from 1 to num_centers. This
                           is that number.
        :type center_num: int

        :param position: (x, y, z) position of the refinement center in the
                         overall grid.
        :type position: list or tuple of float

        :param dx_fine: Resolution at the finest level (near the center).
                        We assume that the three directions have the same
                        resolution.
        :type dx_fine: float

        :param cfl_fine: CFL number at the finest level (near the center).
                           This is dt_finest/dx_finest.
        :type cfl_fine: float

        :param num_levels_with_dt_coarse: How many refinement levels should have
                                     the same dt as the coarsest refinement
                                     level. This will translate to
                                     Carpet::time_refinement_factors. The
                                     minimum is 1 (the coarsest level itself).
        :type num_levels_with_dt_coarse: int

        """
        if not (
            isinstance(refinement_radii, list)
            or isinstance(refinement_radii, tuple)
        ):
            raise TypeError("refinement_radii has to be a list or a tuple")

        if len(refinement_radii) == 0:
            raise ValueError("refinement_radii cannot be empty")

        self.refinement_radii = tuple(sorted(refinement_radii, reverse=True))
        self.center_num = center_num
        self.position = tuple(position)
        self.dx_fine = dx_fine
        self.cfl_fine = cfl_fine
        self.dt_fine = self.cfl_fine * self.dx_fine
        self.num_levels_with_dt_coarse = num_levels_with_dt_coarse

        # If I have only one level evolving with dt_coarse, this has to be the
        # outermost level, so I have no levels synced-up, which means that the
        # CFL factor has to be the same for all the levels, up to the coarsest.
        # Then, for every additional level evolved with dt_coarse, the CFL
        # factor is halved.
        self.cfl_coarse = cfl_fine * 0.5 ** (
            self.num_levels_with_dt_coarse - 1
        )

        self.num_refinement_radii = len(self.refinement_radii)
        self.num_refinement_levels = self.num_refinement_radii + 1

        # If I have 0 levels, then dx_coarse and dx_fine are the same, and then
        # for each level I gain a factor of 2
        self.dx_coarse = self.dx_fine * 2 ** (self.num_refinement_levels - 1)

        self.dt_coarse = self.dx_coarse * self.cfl_coarse

        if not (
            1 <= self.num_levels_with_dt_coarse <= self.num_refinement_levels
        ):
            raise ValueError(
                "num_levels_with_dt_coarse has to be between 1 and "
                f" {self.num_refinement_levels}, "
                f"its value is {self.num_levels_with_dt_coarse}"
            )

    @property
    @lru_cache(1)
    def dx(self):
        """Grid spacing at each refinement level.

        :rtype: tuple of float
        """
        return tuple(
            self.dx_coarse * 0.5**level_num
            for level_num in range(self.num_refinement_levels)
        )

    @property
    @lru_cache(1)
    def dt(self):
        """Timesteps at each refinement level.

        :rtype: tuple of float
        """
        # Here we have to take care of the synced levels.
        #
        # Consider the example with 4 refinement levels, 2 evolving at
        # dt_coarse.
        #
        # Let's set dt_coarse = dt0
        #
        # OB---------RR1----------RR2---------RR3----------RR4---------0
        #      dt0        dt0         dt0/2       dt0/4         dt/8
        #
        # We have a factor 0.5 ** (lvl_num - self.num_levels_with_dt_coarse + 1)
        # because when we have to start dividing by two the first level in which
        # lvl_num = self.num_levels_with_dt_coarse, at which point the exponent
        # has to be 1. We have to start dividing when lvl_num =
        # self.num_levels_with_dt_coarse because we start counting from 0, so in
        # the previous example we had num_levels_with_dt_coarse = 2, so levels 0
        # and 1 evolve with dt0, but level 2 evolves with dt0/2.
        return tuple(
            self.dt_coarse
            if lvl_num < self.num_levels_with_dt_coarse
            else self.dt_coarse
            * 0.5 ** (lvl_num - self.num_levels_with_dt_coarse + 1)
            for lvl_num in range(self.num_refinement_levels)
        )

    @property
    @lru_cache(1)
    def rl_synced_every(self):
        """Number of iterations at which all the refinement levels are at the
        same time.

        :rtype: int

        """
        return 2 ** (
            self.num_refinement_levels - self.num_levels_with_dt_coarse
        )

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        def assign_parameter(param, value, index=None):
            if index is None:
                return f"CarpetRegrid2::{param}_{self.center_num} = {value}"
            return (
                f"CarpetRegrid2::{param}_{self.center_num}[{index}] = {value}"
            )

        ret = []

        ret.append(assign_parameter("num_levels", self.num_refinement_levels))
        ret.append(assign_parameter("position_x", self.position[0]))
        ret.append(assign_parameter("position_y", self.position[1]))
        ret.append(assign_parameter("position_z", self.position[2]))

        for radius_num in range(self.num_refinement_radii):
            ret.append(
                assign_parameter(
                    "radius",
                    self.refinement_radii[radius_num],
                    index=radius_num + 1,  # radii are counted from 1 in Carpet
                )
            )
        return "\n".join(ret)

    def __hash__(self):
        """:py:class:`~.RefinementCenter` is an immutable object."""
        hash_refinement_radii = hash(self.refinement_radii)
        hash_center_num = hash(self.center_num)
        hash_position = hash(self.position)
        hash_dx_fine = hash(self.dx_fine)
        hash_cfl_fine = hash(self.cfl_fine)
        hash_num_levels_with_dt_coarse = hash(self.num_levels_with_dt_coarse)

        # ^ = bitwise xor
        return (
            hash_refinement_radii
            ^ hash_center_num
            ^ hash_position
            ^ hash_dx_fine
            ^ hash_cfl_fine
            ^ hash_num_levels_with_dt_coarse
        )

    @property
    @lru_cache(1)
    def cfl(self):
        """Return dt/dx at every level.

        :returns: Courant factor on every refinement level.
        :rtype: tuple
        """
        return tuple(
            self.dt[ref_level] / self.dx[ref_level]
            for ref_level in range(self.num_refinement_levels)
        )

    @lru_cache(128)
    def dx_at(self, radius: float) -> float:
        """Return dx at the given radius.

        :param radius: Distance from the center.
        :type radius: float

        :returns: Resolution at given point.
        :rtype: float
        """
        # We only care about positive radii
        radius = abs(radius)
        for ref_index, ref_rad in enumerate(self.refinement_radii):
            # If radius == ref_radii, we want the finer dx
            if radius > ref_rad:
                return self.dx[ref_index]
        return self.dx[-1]

    def get_summary(self, outer_boundary=9999.999):
        """Return a summary of all the refinement levels.

        :param outer_boundary: Outer boundary of the grid. Used only for the
                               aesthetics of having a table full of data. If
                               None, a large number full of nines will be used.
        :type outer_boundary: float or None
        """

        refinement_level_boundaries = (
            (outer_boundary,) + self.refinement_radii + (0.0,)
        )

        ret = []

        header = (
            f"Refinement center {self.center_num} "
            f"(x = {self.position[0]: 4.2f},"
            f" y = {self.position[1]: 4.2f}, "
            f" z = {self.position[1]: 4.2f})"
        )
        ret.append(header)

        for ref_level in range(self.num_refinement_levels):
            prev_radius = refinement_level_boundaries[ref_level]
            curr_radius = refinement_level_boundaries[ref_level + 1]
            dx = self.dx[ref_level]
            ret.append(
                f"Ref level {ref_level:2}: "
                f"Radius range = ({prev_radius: 9.3f}, {curr_radius: 9.3f}), "
                f"dx = {dx:.4e} (~ 1/{(1/dx):6.2f}), "
                f"dt = {self.dt[ref_level]:.4e}, "
                f"cfl = {self.cfl[ref_level]:.3f}"
            )
        return "\n".join(ret)


class Grid(BaseThorn):
    """A grid is a collection of :py:class:`~.RefinementCenter` objects. For example,
    a binary black hole simulation will likely have two
    :py:class:`~.RefinementCenter`, centered each black hole.

    This class immutable.

    :ivar outer_boundary: Radius of the outer boundary of the overall grid. We
                          assume that it is the same for the three directions.
    :vartype outer_boundary: float or None

    :ivar num_levels_with_dt_coarse: How many refinement levels should have the
                                     same dt as the coarsest refinement level.
                                     This will translate to
                                     Carpet::time_refinement_factors. The
                                     minimum is 1 (the coarsest level itself).
                                     This is extracted from the refinement center
                                     with the most number of levels.
    :vartype num_levels_with_dt_coarse: int

    :ivar max_num_refinement_levels: Maximum number of refinement levels on the
                                     grid.
    :vartype max_num_refinement_levels: int

    :ivar refinement_centers: Input list, sorted by ``center_num``.
    :vartype refinement_centers: list or tuple of :py:class:`~.RefinementCenter`

    :ivar dx_coarse: Resolution on the coarsest level.
    :vartype dx_coarse: float
    :ivar dt_coarse: Timestep on the coarsest level.
    :vartype dt_coarse: float
    :ivar cfl_coarse: Courant factor on the coarsest level.
    :vartype cfl_coarse: float

    :ivar center_with_most_levels: Refinement center that has the largest number of
                                   levels.
    :vartype center_with_most_levels: :py:class:`~.RefinementCenter`

    :ivar symmetry: If None, no symmetry is not enabled. If 'x', 'y', or 'z',
                    enable reflection symmetry along that axis. If 'xy', 'yz',
                    or 'xz', enable reflection symmetry along two axis. If
                    'xyz', enable octant symmetry. If `90` or `180`, enable
                    rotational symmetry (about the z axis). Rotational symmetry
                    is the kind of symmetry that exists for a single black hole
                    (`90`), or for identical binary (`180`). Activate both
                    rotational and reflection symmetries with `90z` or `180z`.
    :vartype symmetry: str, or None

    :ivar num_ghost: Number of ghost zones.
    :vartype num_ghost: int

    """

    def __init__(
        self,
        refinement_centers,
        outer_boundary,
        num_ghost=3,
        symmetry=None,
        tiny_shift=False,
    ):
        """Constructor.

        The different refinement centers that form a grid must be compatible.
        They must have different ``center_num`` (as they identify different
        refinement centers) and compatible resolutions and timesteps in the
        outer levels. They can have different number of levels.

        :param refinement_centers: List of mesh refinements corresponding to
                                   the various different centers.
        :type refinement_centers: list or tuple of :py:class:`~.RefinementCenter`

        :ivar outer_boundary: Radius of the outer boundary of the overall grid.
                              We assume that it is the same for the three
                              directions.
        :vartype outer_boundary: float or None

        :param symmetry: If None, no symmetry is not enabled. If 'x', 'y', or
                             'z', enable reflection symmetry along that axis. If
                             'xy', 'yz', or 'xz', enable reflection symmetry
                             along two axis. If 'xyz', enable octant symmetry.
                             If `90` or `180`, enable rotational symmetry (about
                             the z axis). Rotational symmetry is the kind of
                             symmetry that exists for a single black hole
                             (`90`), or for identical binary (`180`). Activate
                             both rotational and reflection symmetries with
                             `90z` or `180z`.
        :vartype symmetry: str, or None

        :param tiny_shift: Apply a tiny (subpixel) shift to the outer boundary so
                           that the point (0,0,0) is not on the grid. This is
                           hard-coded to be 1/7 of a pixel.
        :type tiny_shift: bool

        :param num_ghost: Number of ghost zones.
        :type num_ghost: int

        """
        if not all(
            (isinstance(m, RefinementCenter) for m in refinement_centers)
        ):
            raise TypeError("Input is not a list of RefinementCenter")

        self.num_centers = len(refinement_centers)

        # Check that all the grids have different center_num
        centers = {ref_center.center_num for ref_center in refinement_centers}
        if len(centers) != self.num_centers:
            raise ValueError("Grids do not have all different IDs")

        # The grids are compatible if the have the same dx_coarse, dt_coarse,
        # and num_levels_with_dt_coarse

        dx_coarse = {ref_center.dx_coarse for ref_center in refinement_centers}
        if len(dx_coarse) != 1:
            raise ValueError(
                "Grids have different resolution at the coarsest level"
            )

        dt_coarse = {ref_center.dt_coarse for ref_center in refinement_centers}
        if len(dt_coarse) != 1:
            raise ValueError(
                "Grids have different timesteps at the coarsest level"
            )

        # Tuple unpacking on a set with one single element extract the only
        # element
        (self.dx_coarse,) = dx_coarse
        (self.dt_coarse,) = dt_coarse

        self.cfl_coarse = self.dt_coarse / self.dx_coarse

        self.max_num_refinement_levels = max(
            ref_center.num_refinement_levels
            for ref_center in refinement_centers
        )

        # Now we go level by level and check that they all have the same timestep
        for level in range(self.max_num_refinement_levels):
            dts = {
                ref_center.dt[level]
                for ref_center in refinement_centers
                if level < ref_center.num_refinement_levels
            }
            if len(dts) != 1:
                raise ValueError(
                    "Grids have different timesteps at level {level}"
                )

        # num_levels_with_dt_coarse is the one of the center with the most levels
        self.center_with_most_levels = sorted(
            refinement_centers, key=lambda m: m.num_refinement_levels
        )[-1]

        self.num_levels_with_dt_coarse = (
            self.center_with_most_levels.num_levels_with_dt_coarse
        )

        self.refinement_centers = tuple(
            sorted(refinement_centers, key=lambda m: m.center_num)
        )

        self.outer_boundary = outer_boundary
        self.tiny_shift = tiny_shift

        _allowed_symmetries = (
            "x",
            "y",
            "z",
            "xy",
            "yz",
            "xz",
            "xyz",
            "90",
            "90z",
            "180",
            "180z",
            None,
        )

        if symmetry not in _allowed_symmetries:
            raise ValueError(
                f"symmetry has to be one in {_allowed_symmetries}"
            )

        self.symmetry = symmetry

        self.num_ghost = num_ghost

    @property
    @lru_cache(1)
    def time_refinement_factors(self):
        """Return parameters for Carpet to set up time_refinement_factors."""

        time_refinement = []
        time_refinement.append('"[')

        # Let us consider an example with four levels (two radii), two levels
        # synced. In this case, max_levels = 4, synced_levels = 2. What we want
        # as output is "[1,1,2,4]". So, we have to append '1' twice, one for
        # each synced_levels, then 2 and 3. What we do is we define a variable
        # `power` that takes the value 0 when we are in the "synced_levels"
        # region, and then it grows after that.

        power = 0

        for level in range(self.max_num_refinement_levels):
            if level >= self.num_levels_with_dt_coarse:
                power += 1
            time_refinement.append(f"{2**power}")
            time_refinement.append(",")

        # Adding a comma at each iteration leaves us with a trailing comma that
        # we don't want. Instead, we want the closing bracket.
        time_refinement[-1] = ']"'

        return "".join(time_refinement)

    @property
    @lru_cache(1)
    def rl_synced_every(self):
        """Number of iterations at which all the refinement levels are at the
        same time.

        :rtype: int

        """
        return self.center_with_most_levels.rl_synced_every

    @property
    @lru_cache(1)
    def parfile_code(self):
        """Return the code you would put in your parfile."""

        def assign_parameter(param, value, thorn="Carpet"):
            return f"{thorn}::{param} = {value}"

        if self.tiny_shift:
            dx_fine = self.dx_coarse / 2.0 ** (
                self.max_num_refinement_levels - 1
            )
            # We add a 1/7 of a grid point
            outer_boundary_plus = self.outer_boundary + dx_fine / 7.0
            outer_boundary_minus = -self.outer_boundary + dx_fine / 7.0
        else:
            outer_boundary_plus = self.outer_boundary
            outer_boundary_minus = -self.outer_boundary

        ret = []

        ret.append(
            f"""\
CartGrid3D::type = "coordbase"
Carpet::domain_from_coordbase = "yes"

Driver::ghost_size = {self.num_ghost}
CoordBase::boundary_size_x_lower = {self.num_ghost}
CoordBase::boundary_size_y_lower = {self.num_ghost}
CoordBase::boundary_size_z_lower = {self.num_ghost}
CoordBase::boundary_size_x_upper = {self.num_ghost}
CoordBase::boundary_size_y_upper = {self.num_ghost}
CoordBase::boundary_size_z_upper = {self.num_ghost}

ReflectionSymmetry::avoid_origin_x = no
ReflectionSymmetry::avoid_origin_y = no
ReflectionSymmetry::avoid_origin_z = no

CoordBase::domainsize = "minmax"
CoordBase::xmax = {outer_boundary_plus}
CoordBase::ymax = {outer_boundary_plus}
CoordBase::zmax = {outer_boundary_plus}
CoordBase::dx = {self.dx_coarse}
CoordBase::dy = {self.dx_coarse}
CoordBase::dz = {self.dx_coarse}"""
        )

        # When it comes to the lower boundaries, some symmetries are equivalent:
        # 90 ~ xy, 90z ~ xyz, 180 ~ x, 180z ~ xz

        _map = {"90": "xy", "90z": "xyz", "180": "x", "180z": "xz", None: ""}

        # What is this funny business?

        # _map translates the rotational symmetries to reflection ones (in terms
        # of boundaries). Then, we set effective_symmetry_for_boundary to
        # something that has always the format of a reflection symmetry (ie,
        # some letters). This is done by either getting the value out of the
        # _map dictionary (when the symmetry is a rotational one, so it is in
        # _map), or the symmetry itself (when the symmetry is already reflection)
        effective_symmetry_for_boundary = _map.get(
            self.symmetry, self.symmetry
        )

        for dim in ("x", "y", "z"):
            if dim in effective_symmetry_for_boundary:
                ret.append(assign_parameter(f"{dim}min", 0, thorn="CoordBase"))
                ret.append(
                    assign_parameter(
                        f"boundary_shiftout_{dim}_lower", 1, thorn="CoordBase"
                    )
                )
            else:
                ret.append(
                    assign_parameter(
                        f"{dim}min", outer_boundary_minus, thorn="CoordBase"
                    )
                )

        # Will add a blank line
        ret.append("")

        if self.symmetry in ("90", "90z"):
            ret.append('ActiveThorns = "RotatingSymmetry90"')
            ret.append("CarpetRegrid2::symmetry_rotating90 = yes\n")
        elif self.symmetry in ("180", "180z"):
            ret.append('ActiveThorns = "RotatingSymmetry180"')
            ret.append("CarpetRegrid2::symmetry_rotating180 = yes\n")

        # Add reflection
        if self.symmetry is not None:
            for dim in ("x", "y", "z"):
                if dim in self.symmetry:
                    ret.append(
                        f"ReflectionSymmetry::reflection_{dim}     = yes"
                    )
            # Will add a blank line
            ret.append("")

        ret.append(
            assign_parameter(
                "num_centres", self.num_centers, thorn="CarpetRegrid2"
            )
        )

        ret.append(
            assign_parameter(
                "max_refinement_levels", self.max_num_refinement_levels
            )
        )
        ret.append(
            assign_parameter(
                "time_refinement_factors", self.time_refinement_factors
            )
        )

        # Will add a blank line
        ret.append("")

        for center in self.refinement_centers:
            ret.append(center.parfile_code)
        return "\n".join(ret)
