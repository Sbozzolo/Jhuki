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

The available functions are:

- :py:func:`~.simfactory_option`, which helps setting up options that are
controllable by SimFactory.

- :py:func:`~.write_parfile`, which write the `.par` file from a template
when the file is executed as `.rpar` by SimFactory.

"""
import inspect  # Needed to capture locals() of the caller to write the parfile
import re
import sys

from jhuki.parfile import write_one_parfile_from_template


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

    We convert all the numbers to floats, so if you need integers, you have to
    explicitly cast the variable to int, for example:

    .. code-block::python

        num_ref_radii = int(simfactory_option("@REFLEVELS@", 10))

    :param name: Name of the SimFactory define key. It has to be with surrounded
                 by @, for example, @XMAX@.
    :type name: str
    :param default: Default value that the variable has to have if define is not
                    specified.
    :type default: int, float, or str

    :returns: The default value if the define key is not passed to SimFactory,
              otherwise the value defined with SimFactory.
    :rtype: float or str (integers are converted)

    """
    if not isinstance(name, str):
        raise TypeError(f"name has to a string, it is {name}")

    # If default is a string, we will keep it a string. If it is an int, we
    # convert it to float. Often one sets parameter like par_b = 6, but here 6
    # has to be intended as a float, or it can cause troubles for other
    # calculations.
    type_ = float if isinstance(default, int) else type(default)

    return type_(name) if name[0] != "@" else type_(default)


def write_parfile(template: str, file_name=None) -> None:
    """Write parfile from given template using locals() as substitution
    dictionary.

    The variables in template get substituted with the values in the local scope
    of the calling function. So, you can safely put this function at the end of
    your rpar file to have the template substituted with the variables you have
    defined there.

    :param template: Template of parfile to write. See
                     :py:func:`~.write_one_parfile_from_template`
    :type template: str
    :param file_name: Name of the file that is calling this function.
                      (Used only for debug purposes.)
    :type file_name: str

    """
    # Compact notation to say: if file_name is not None, file_name = file_name,
    # otherwise it is sys.argv[0]
    file_name = file_name or sys.argv[0]

    # The regular expression matches something that has the form FILENAME.rpar
    rpar_rx = re.compile(r"^(.*)\.rpar$")

    if not rpar_rx.match(file_name):
        raise RuntimeError("Function was not called by a rpar file")

    # Replace "file_name.rpar" with "file_name.par"
    output_file_name = rpar_rx.sub(r"\1.par", file_name)

    # We want to capture the locals() of the caller (to use it as substitution
    # dictionary). This function should be called directly from the `rpar` file

    # From https://stackoverflow.com/a/62885157
    frame = inspect.currentframe().f_back
    try:
        # frame.f_locals is the locals() of the calling function
        subs = frame.f_locals

        # From here, we have to remove the template variable itself (otherwise
        # we would have a template-ception)
        #
        # The easiest way to do this is to match the value in the dictionary
        # subs. This is not particularly efficient, but we expect subs to be
        # small
        subs_no_template = {k: v for k, v in subs.items() if v != template}

        write_one_parfile_from_template(
            template, sub_dict=subs_no_template, out_file=output_file_name
        )
    finally:
        # Ensure that we release the resources
        del frame
