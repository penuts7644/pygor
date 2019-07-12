# Create IGoR models and calculate the generation probability of V(D)J and
# CDR3 sequences. Copyright (C) 2019 Wout van Helvoirt

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


"""Contains IgorInterface class for interfacing with IGoR's commandline tool."""


import shlex
from subprocess import Popen, PIPE


class IgorInterface(object):
    """Executes IGoR commands via new commandline subprocess.

    Parameters
    ----------
    command : list
        A list with strings and nested lists that will be build into a
        subprocess command.

    Methods
    -------
    call()
        Call the IGoR program and return stdout and stderr messages.
    get_command()
        Returns the created commandline subprocess string.
    set_command(args)
        Set a new commandline string using a nested list.

    """
    def __init__(self, command):
        super(IgorInterface, self).__init__()
        self.command = self._subprocess_builder(options=command)

    def _subprocess_builder(self, options, level=0):
        """Creates a subprocess command string from an sorted input list.

        Parameters
        ----------
        options : list
            A Python nested sorted list with each value being a
            command/options/flag. Items within a list are separated by a white
            space in the output string. The depth of the nested lists will
            determine the number '-' characters to add in front of the first
            element for each list.
        level : int, optional
            The initial start depth level indication the number of '-'
            characters to append to the first item in the lists. (default: 0)

        Returns
        -------
        str
            The formatted commandline subprocess as string.

        """
        # Create the commandline string from the options in the list.
        command_str = ""
        for index, val in enumerate(options):
            if isinstance(val, list):
                val = self._subprocess_builder(val, level=level+1)
            if isinstance(val, str):
                if index == 0:
                    val = (level * '-') + val
                command_str += ' ' + val
        return command_str.strip(' ')

    def call(self):
        """Calls IGoR via commandline via the subprocess command string.

        Returns
        -------
        tuple
            A tuple containing the exit code, standard out, standard error and
            the executed command as string.

        """
        # Execute the commandline process and return the results.
        updated_command = 'igor ' + self.command
        process = Popen(shlex.split(updated_command),
                        stderr=PIPE, stdout=PIPE)
        (stdout, stderr) = process.communicate()
        return (process.returncode, stdout, stderr, updated_command)

    def get_command(self):
        """Collect and returns the string formatted IGoR command.

        Returns
        -------
        str
            A command formatted as string to be used by the call function.

        """
        return self.command

    def set_command(self, command):
        """Sets a new IGoR command by processing teh input lists into a string.

        Parameters
        ----------
        command : list
            A list with strings and nested lists that will be build into a
            subprocess command.

        """
        self.command = self._subprocess_builder(options=command)
