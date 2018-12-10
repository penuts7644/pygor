# ImmunoProbs Python package uses a simplified manner for calculating the
# generation probability of V(D)J and CDR3 sequences.
# Copyright (C) 2018 Wout van Helvoirt

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""Executable file for running functions located in immuno_probs.cli directory."""


import argparse

from immuno_probs.cli.create_cdr3_anchors import CreateCdr3Anchors
from immuno_probs.util.cli import dynamic_cli_options
from immuno_probs.util.constant import set_max_threads
from immuno_probs.util.constant import set_separator


def main():
    """Function to create the ArgumentParser containing the sub-options."""
    # Create the parser with general commands and set the subparser.
    description = 'ImmunoProbs Python package uses a simplified manner for calculating ' \
        'the generation probability of V(D)J and CDR3 sequences.'
    parser_general_options = {
        '--separator': {
            'type': 'str',
            'nargs': '?',
            'help': 'The separator character to be used when writing files ' \
                    '(default: comma character).'
        },
        '--max-threads': {
            'type': 'int',
            'nargs': '?',
            'help': 'The max number of threads the program is allowed to use ' \
                    '(default: max threads).'
        }
    }
    parser = argparse.ArgumentParser(prog='immuno-probs', description=description)
    parser = dynamic_cli_options(parser=parser, options=parser_general_options)
    subparsers = parser.add_subparsers(help='Supported immuno-probs options, command plus help ' \
                                            'displays more information for the option.',
                                       dest='subparser_name')

    # Add main- and suboptions to the subparser.
    cca = CreateCdr3Anchors(subparsers=subparsers)

    # Parse the commandline arguments, set variables and execute correct function.
    parsed_arguments = parser.parse_args()
    if parsed_arguments.separator is not None:
        set_separator(parsed_arguments.separator)
    if parsed_arguments.max_threads is not None:
        set_max_threads(parsed_arguments.max_threads)

    if parsed_arguments.subparser_name == cca.__class__.__name__:
        cca.run(args=parsed_arguments)
    else:
        print("No option specified, run 'immuno-probs -h' to display all options.")


if __name__ == '__main__':
    main()
