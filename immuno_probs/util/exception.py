# ImmunoProbs Python package uses simplified manner for calculating the
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


"""Contains custom exception classes used in immuno_probs."""


class CharacterNotFoundException(Exception):
    """Exception for when a character has not been found."""
    def __init__(self, message, character):
        super(CharacterNotFoundException, self).__init__()
        self.message = message
        self.character = character

    def __str__(self):
        return "{0}: '{1}'".format(self.message, self.character)


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


class MaxThreadsValueException(Exception):
    """Exception when the MAX_THREADS variable is faulty."""
    def __init__(self, message, character):
        super(MaxThreadsValueException, self).__init__()
        self.message = message
        self.character = character

    def __str__(self):
        return "{0}: '{1}'".format(self.message, self.character)

class SeparatorNotValidException(Exception):
    """Exception when the SEPARATOR variable is faulty."""
    def __init__(self, message, character):
        super(SeparatorNotValidException, self).__init__()
        self.message = message
        self.character = character

    def __str__(self):
        return "{0}: '{1}'".format(self.message, self.character)
