[![GPLv3
license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

# PreCactus

Writing parameter files for Cactus simulations is not easy. To achieve a
successful evolution, one has to tune several parameters, which typically depend
on the grid configuration or on other settings. Writing these par files by hand
is tedious and error-prone. `PreCactus` is a Python library to prepare
simulations for the `Einstein Toolkit` (or Cactus-based codes).

## Features

* Generate parameter files from a template using configuration files or
  command-line arguments
* Take care of the grid configuration given the desired resolution at the finest
  level and other details


## Examples

### Working with grid structures

Problem: you want to generate the parfile code for a grid structure with two
refinement centers, each with 8 levels and resolution at the finest level 0.001
and CFL factor of 0.4 (in the finest level). In this, you want to ensure that
the maximum timestep on the grid never exceeds 1 to avoid numerical instability.

``` python
#!/usr/bin/env python3

from precactus import grid as pg

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

## Installation

The best way to install `PreCactus` is by cloning this repo and using
[poetry](https://python-poetry.org/). If you have poetry install, just run
`poetry install` in the folder where you cloned the repo to install `PreCactus`.

## Tests

`PreCactus` comes with a suite of unit tests. To run the tests,
```sh
poetry run pytest --cov=./ --cov-report=term
```
Tests are automatically run after each commit by GitHub Actions. This will also
tell you what is the test coverage.

# Code style

- We lint the code with `black -l 79`.
