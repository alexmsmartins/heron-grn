"""
prol-evolution GRN, a model for regulatory gene networks

Copyright (C) 2008 Luis Pureza, Oseias Santos, Pedro Martins and
Ricardo Pereira.
 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


class GRNElement(object):
    """
    Base class for all GRN graph elements
    """

    count = 0

    def __init__(self, parent, id=None):
        """
        Creates a new GRN element with the given parent

        Arguments:
        - `parent`: The parent node of this element
        """
        self.parent = parent
        self.__class__.count += 1
        self.id = id
        if not self.id:
            self.id = "%s_%d" % (type(self).__name__, self.__class__.count)
        self.enabled = False

    def __repr__(self):
        return "GRN(%s)" % self.id

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id.__hash__()
        
