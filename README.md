<p align="center">
<img src="https://github.com/Sbozzolo/jhuki/raw/master/logo.png" height="250">
</p>

[![GPLv3
license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)
![Tests](https://github.com/Sbozzolo/Jhuki/workflows/Tests/badge.svg)


Writing parameter files for Cactus simulations is not easy. To achieve a
successful evolution, one has to tune several parameters, which typically depend
on the grid configuration or on other settings. Writing these par files by hand
is tedious and error-prone. `Jhuki` is a Python library to prepare
simulations for the `Einstein Toolkit` (or Cactus-based codes).

> :warning: This package is currently under development. It may be full of bugs,
>           and its interfaces might change without notice. It is also highly
>           opinionated.

## Features

* Generate parameter files from a template using configuration files or
  command-line arguments
* Take care of the grid configuration given the desired resolution at the finest
  level and other details
* Generate binary black hole configurations for quasi-circular mergers
  (automatically setting the linear momenta)

## Examples

### Working with grid structures

Problem: you want to generate the parfile code for a grid structure with two
refinement centers, each with 8 levels and resolution at the finest level 0.001
and CFL factor of 0.4 (in the finest level). In this, you want to ensure that
the maximum timestep on the grid never exceeds 1 to avoid numerical instability.

``` python
#!/usr/bin/env python3

from jhuki import grid as pg

refinement_radii = tuple(2**level for level in range(7))

center1 = pg.RefinementCenter(refinement_radii,
                              dx_fine=0.001,
                              cfl_fine=0.5,
                              center_num=1,
                              position=(10,0,0))

# Same but with different center_num and position
center2 = pg.RefinementCenter(refinement_radii,
                              dx_fine=0.001,
                              cfl_fine=0.5,
                              center_num=2,
                              position=(-10,0,0))

grid_not_synced = pg.Grid((center1, center2), outer_boundary=1000)
grid_synced = pg.set_dt_max_grid(grid_not_synced, dt_max=1)

print(grid_synced.parfile_code)
```
This will output

``` sh
CartGrid3D::type = "coordbase"
Carpet::domain_from_coordbase = "yes"
CoordBase::domainsize = "minmax"
CoordBase::xmin = 1000
CoordBase::ymin = 1000
CoordBase::zmin = 1000
CoordBase::xmax = 1000
CoordBase::ymax = 1000
CoordBase::zmax = 1000
CoordBase::dx = 0.64
CoordBase::dy = 0.64
CoordBase::dz = 0.64
Carpet::max_refinement_levels = 8
Carpet::time_refinement_factors = "[1,1,2,4,8,16,32,64]"
CarpetRegrid2::num_levels_1 = 8
CarpetRegrid2::position_x_1 = 10
CarpetRegrid2::position_y_1 = 0
CarpetRegrid2::position_z_1 = 0
CarpetRegrid2::radius_1[1] = 64
CarpetRegrid2::radius_1[2] = 32
CarpetRegrid2::radius_1[3] = 16
CarpetRegrid2::radius_1[4] = 8
CarpetRegrid2::radius_1[5] = 4
CarpetRegrid2::radius_1[6] = 2
CarpetRegrid2::radius_1[7] = 1
CarpetRegrid2::num_levels_2 = 8
CarpetRegrid2::position_x_2 = -10
CarpetRegrid2::position_y_2 = 0
CarpetRegrid2::position_z_2 = 0
CarpetRegrid2::radius_2[1] = 64
CarpetRegrid2::radius_2[2] = 32
CarpetRegrid2::radius_2[3] = 16
CarpetRegrid2::radius_2[4] = 8
CarpetRegrid2::radius_2[5] = 4
CarpetRegrid2::radius_2[6] = 2
CarpetRegrid2::radius_2[7] = 1
```

You can also add a small shift to the grid so that the origin is not on (0,0,0)
passing the `tiny_shift` argument to `Grid`.

Another useful method is `rl_synced_every`, which returns the number of
iterations at which all the refinement levels are at the same time.


### Working with binary black holes

Problem: you want to simulate a binary black quasi-circular merger. This can be
tricky because you have to provide the correct linear momenta. With `Juhki`,
this is trivial:

``` python
#!/usr/bin/env python3

from jhuki.twopunctures import prepare_quasicircular_inspiral

mass_ratio = 1
coordinate_distance = 12
dimensionless_spin_plus = (0.1, 0.2, 0.3)
dimensionless_spin_minus = (0.4, -0.1, -0.2)

twopunctures = prepare_quasicircular_inspiral(mass_ratio,
                                              coordinate_distance,
                                              dimensionless_spin_plus,
                                              dimensionless_spin_minus)
print(twopunctures.parfile_code)
```

This will an output similar to the following:

```
ADMBase::initial_data = twopunctures
ADMBase::initial_lapse = twopunctures-averaged
ADMBase::initial_shift = zero
ADMBase::initial_dtlapse = zero
ADMBase::initial_dtshift = zero
TwoPunctures::give_bare_mass = no
TwoPunctures::par_b = 6.0
TwoPunctures::target_m_plus = 0.5
TwoPunctures::target_m_minus = 0.5
TwoPunctures::par_P_plus[0] = -0.000531433937072403
TwoPunctures::par_P_plus[1] = 0.0848055056299618
TwoPunctures::par_P_plus[2] = 0
TwoPunctures::par_P_minus[0] = 0.000531433937072403
TwoPunctures::par_P_minus[1] = -0.0848055056299618
TwoPunctures::par_P_minus[2] = 0
TwoPunctures::par_S_plus[0] = 0.025
TwoPunctures::par_S_plus[1] = 0.05
TwoPunctures::par_S_plus[2] = 0.075
TwoPunctures::par_S_minus[0] = 0.1
TwoPunctures::par_S_minus[1] = -0.025
TwoPunctures::par_S_minus[2] = -0.05
TwoPunctures::center_offset[0] = 0.0
TwoPunctures::center_offset[1] = 0.0
TwoPunctures::center_offset[2] = 0.0
```

This module should be used along with the `Grid` one.


## Installation

The best way to install `Jhuki` is by cloning this repo and using
[poetry](https://python-poetry.org/). If you have poetry install, just run
`poetry install` in the folder where you cloned the repo to install `Jhuki`.
Alternatively, `Jhuki` is available on PyPI.

## Tests

`Jhuki` comes with a suite of unit tests. To run the tests,
```sh
poetry run pytest --cov=./ --cov-report=term
```
Tests are automatically run after each commit by GitHub Actions. This will also
tell you what is the test coverage.

# Code style

- We lint the code with `black -l 79`.

# What does _jhuki_ mean?

The word _jhuki_ belongs to the Tohono O'odham vocabulary and means *rain*. If
[kuibit](https://githum.com/Sbozzolo/kuibit) is the tool you use to collect the
fruit of your `Cactus` simulations, then `jhuki` is what allowed that fruit to
be there in the first place.

## Credits

The logo contains elements designed by [pngtree.com](pngtree.com).

The computation of the momenta for quasi-circular mergers of binary black holes
uses
[NRPyPN](https://einsteintoolkit.org/thornguide/EinsteinInitialData/NRPyPN/documentation.html).
If you use this module, please follow the citation guidelines as specified by
the documentation in the `NRPyPN` repo.
