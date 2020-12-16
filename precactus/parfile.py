#!/usr/bin/env python3

# Copyright (C) 2020 Gabriele Bozzola
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


"""The :py:mod:`~.parfile` module provides functions to programmatically write
parameter files.
"""

from string import Template


def write_one_parfile_from_template(template, sub_dict, out_file):
    """Write a parfile to out_file from template with substitution
    specified by the dictionary sub_dict.

    We add an header with the list of variables that were introduced.

    :param template: Basic template with parameters that have to be
                     substituted. It follows Python's template
                     conventions, so the variables to be substituted
                     are to be prefixed with a dollar sign.
    :type template: str
    :param sub_dict: Dictionary that maps entries that have to be
                     substituted with their value.
    :type sub_dict: dict
    :param out_file: Path of the output file
    :type out_file: str
    """

    # We add an header to save what variables were substituted
    substituted_variables = []
    for key, val in sorted(sub_dict.items()):
        # We cast everything to string and substitute all the newlines
        # with commented newlines so that we can ensure that multiline
        # parameters are properly commented.
        val = str(val).replace("\n", "\n# ")
        substituted_variables.append(f"# {key} = {val}")

    substituted_variables_str = "\n".join(substituted_variables)

    header = f"""\
################################################################################
# Variables introduced by PreCactus
################################################################################

{substituted_variables_str}
"""

    template_with_header = header + template

    with open(out_file, "w") as file_:
        file_.write(Template(template_with_header).substitute(sub_dict))
