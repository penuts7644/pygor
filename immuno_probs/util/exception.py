# ImmunoProbs Python package able to calculate the generation probability of
# V(D)J and CDR3 sequences. Copyright (C) 2018 Wout van Helvoirt

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


"""Contains custom exception classes used in immuno_probs."""


class AlignerException(Exception):
    """Exception for when an error occured with one of the aligner programs."""
    def __init__(self, aligner_message):
        super(AlignerException, self).__init__()
        self.message = aligner_message

    def __str__(self):
        return "An error occured with the aligner:\n{0}".format(self.message)


class CharacterNotFoundException(Exception):
    """Exception for when a character has not been found."""
    def __init__(self, message, value):
        super(CharacterNotFoundException, self).__init__()
        self.message = message
        self.value = value

    def __str__(self):
        return "{0}: '{1}'".format(self.message, self.value)


class DirectoryNonExistingException(Exception):
    """Exception when the WORKING_DIR variable is faulty."""
    def __init__(self, message, value):
        super(DirectoryNonExistingException, self).__init__()
        self.message = message
        self.value = value

    def __str__(self):
        return "{0}: '{1}'".format(self.message, self.value)


class GeneIdentifierException(Exception):
    """Exception for incorrect gene identifier value."""
    def __init__(self, message, identifier):
        super(GeneIdentifierException, self).__init__()
        self.message = message
        self.identifier = identifier

    def __str__(self):
        return "{0}: '{1}''".format(self.message, self.identifier)


class IndexNotFoundException(Exception):
    """Exception for being unable to collect and index from a list."""
    def __init__(self, message, identifier):
        super(IndexNotFoundException, self).__init__()
        self.message = message
        self.identifier = identifier

    def __str__(self):
        return "{0}: '{1}''".format(self.message, self.identifier)


class NumThreadsValueException(Exception):
    """Exception when the NUM_THREADS variable is faulty."""
    def __init__(self, message, value):
        super(NumThreadsValueException, self).__init__()
        self.message = message
        self.value = value

    def __str__(self):
        return "{0}: '{1}'".format(self.message, self.value)

class SeparatorNotValidException(Exception):
    """Exception when the SEPARATOR variable is faulty."""
    def __init__(self, message, value):
        super(SeparatorNotValidException, self).__init__()
        self.message = message
        self.value = value

    def __str__(self):
        return "{0}: '{1}'".format(self.message, self.value)


class SubprocessException(Exception):
    """Exception for when an error occured with a subprocess execution."""
    def __init__(self, subprocess_message):
        super(SubprocessException, self).__init__()
        self.message = subprocess_message

    def __str__(self):
        return "An error occured with the subprocess:\n{0}".format(self.message)
