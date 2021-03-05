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


"""The :py:mod:`~.parfile` module provides functions to programmatically write
parameter files. The main functions provided are:

- :py:func:`~.write_one_parfile_from_template`. This function takes a Python
  string template, a dictionary with the instructions to perform the
  substitution, and the name of the file where to write the output.
- :py:func:`~.write_many_parfiles_from_template`. This function calls the
  previous function (:py:func:`~.write_one_parfile_from_template`) multiple times
  on a multiple substitution lists to produce multiple files.

"""

import logging
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
# Variables introduced by Jhuki
################################################################################

{substituted_variables_str}

"""

    template_with_header = header + template

    with open(out_file, "w") as file_:
        file_.write(Template(template_with_header).substitute(sub_dict))


def write_many_parfiles_from_template(template, subs_dict, out_file_prefix):
    """Write several parfiles from the template with substitution specified by the
    dictionary sub_dict_list. The dictionary has to have as keys the variables to
    be substituted and as values:
    - "scalars" for those values that have to be the same in all the parfiles
    - lists, with one entry per file, for those that have to be change file by file.

    For example: {'dx': [1, 2], 'dy': 3} will create two files with
    'dx = 1' and 'dx = 2' and 'dy = 3' in both cases.

    We add an header with the list of variables that were introduced.

    :param template: Basic template with parameters that have to be
                     substituted. It follows Python's template
                     conventions, so the variables to be substituted
                     are to be prefixed with a dollar sign.
    :type template: str
    :param subs_dict: Dictionary that maps entries that have to be
                     substituted with their value, one per each file.
                     The values of the dictionary have to be lists or "scalars"
                     (if they have to have the same value for the files).
    :type subs_dict: dict with lists as values
    :param out_file_prefix: Path of the output files, without the extension ".par".
                            The parfiles will be saved as ``{out_file_prefix}1.par``,
                            ``{out_file_prefix}2.par``, and so on.
    :type out_file: str

    """

    # We check that all the lists with more than one element have the same length
    list_lengths = [
        len(v)
        for v in subs_dict.values()
        if (isinstance(v, list) and len(v) > 1)
    ]

    # Here we have greater than one because we allow for 0 elements
    if len(set(list_lengths)) > 1:
        raise ValueError("Substitution lists have different lengths.")

    if list_lengths:
        number_of_parfiles = list_lengths[0]
    else:
        number_of_parfiles = 1
    logging.debug(f"Writing {number_of_parfiles} parfiles")

    # Next, we ensure that all the elements are lists
    sub_dict_list = {
        k: (v if isinstance(v, list) else [v]) for k, v in subs_dict.items()
    }

    # Now we "broadcast" the scalar elements to list. With this,
    # sub_dict_broadcast will have lists with the same length
    # (= the number of parfiles)
    sub_dict_broadcast = {
        # else here is when the lists have only one element
        k: (v if len(v) == number_of_parfiles else v * number_of_parfiles)
        for k, v in sub_dict_list.items()
    }

    # Now we loop over the sub_dict_broadcast and create one parfile
    # for each entry
    for num_parfile in range(number_of_parfiles):
        file_name = f"{out_file_prefix}{num_parfile}.par"
        logging.debug(f"Writing {file_name}")
        subs = {k: v[num_parfile] for k, v in sub_dict_broadcast.items()}
        write_one_parfile_from_template(template, subs, file_name)
