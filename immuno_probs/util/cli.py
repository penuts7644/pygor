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


"""Contains command-line functions used in immuno_probs."""


def dynamic_cli_options(parser, options):
    """Semi-dynamically adds options to the given parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        ArgumentParser to append options to.
    options : dict
        A Python dict with key being the full name of the option. The value is
        a dict that corresponds to input arguments of the
        ArgumentParser.add_argument function. Note: type argument values must be
        surrounded by quotes.

    Returns
    -------
    ArgumentParser
        Object containing the expected commandline arguments. Still needs to
        parse the commandline arguments.

    """
    # Semi-dynamically create the argparse arguments from given inputs.
    for name, kwargs in options.iteritems():
        kwargs_str = ""
        for (option, value) in kwargs.iteritems():
            if isinstance(value, str) and not option == 'type':
                kwargs_str += ', {}="{}"'.format(option, value)
            else:
                kwargs_str += ', {}={}'.format(option, value)
        eval('parser.add_argument("{0}"{1})'.format(name, kwargs_str))

    # Return the updated parser.
    return parser


def subprocess_builder(options, level=0):
    """Creates a subprocess command-line string from ordered input dictionary.

    Parameters
    ----------
    options : list
        A Python nested ordered list with each value being a command/option.
        White spaces will seperate the options. Note: the depth of the nested
        lists will determine the number '-' characters to add in front of the
        first element for each list.
    level : int, optional
        The initial start depth level for indication the number of '-' characters
        to add to the first item in the lists. (default: 0)

    Returns
    -------
    str
        The formatted command-line command as string.

    """
    # Create the commandline string from the options in the list.
    command_str = ""
    for index, val in enumerate(options):
        if isinstance(val, list):
            val = subprocess_builder(val, level=level+1)
        if isinstance(val, str):
            if index == 0:
                val = (level * '-') + val
            command_str += ' ' + val
    return command_str.strip(' ')
