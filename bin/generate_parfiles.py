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


"""This script takes a template of a parfile and substitution list from command
line and produces new parfiles. The template has to follow the same rules as
Python Templates: the strings that have to be substituted are prefixed with the
dollar sign.

"""

import os.path
from argparse import ArgumentError

import configargparse

from jhuki import parfile

if __name__ == "__main__":

    parser = configargparse.ArgParser()
    parser.add(
        "-c", "--config", is_config_file=True, help="Path of the config file"
    )
    parser.add(
        "-t",
        "--template",
        required=True,
        help="Path of the template",
    )

    parser.add(
        "-o",
        "--output",
        default="output",
        help="Path of where to save the files including the prefix of the parfile name",
    )

    # HACK: Here we add all the unknown arguments to the parser,
    #       so that can parse them later
    parsed, unknown = parser.parse_known_args()
    for arg in unknown:
        # All the arguments defined in configuration files are pre-pended with
        # '--' and the value is after a '='
        if arg.startswith("--"):
            arg_name, _ = arg.split("=")
            try:
                # We need nargs='*' to support lists.
                # Using nargs='*' we turn everything into list.
                parser.add(arg_name, nargs="*")
            except ArgumentError:
                pass

    args = parser.parse_args()

    # We don't want in our substitution list some keys that are used for other
    # stuff, like reading the configuration file
    substitution_list = {
        k: v
        for k, v in vars(args).items()
        if k not in ("config", "template", "output")
    }

    if not (os.path.exists(args.template) and os.path.isfile(args.template)):
        raise ValueError(f"Invalid file {args.template}")

    # Read file as string
    with open(args.template, "r") as file_:
        template_str = file_.read()

    # Finally, write the parfiles
    parfile.write_many_parfiles_from_template(
        template_str, substitution_list, args.output
    )
