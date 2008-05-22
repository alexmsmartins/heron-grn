"""
prol-evolution GRN, a model for regulatory gene networks

Copyright (C) 2008 Luis Pureza, Oseias Santos, Pedro Matrins and
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


import re
from grn_element import *
from defs import *
from protein import *
from util import *

class MessengerRNA (GRNElement):
    """ Messenger RNA """

    def __init__ (self, parent, sequence):
        """ Creates new messenger RNA """
	GRNElement.__init__(self, parent)
        self.sequence = sequence

    def translate (self):
        """ Translate the mRNA into proteins """

        # Find the start codon
        index_of_MET = translation_table.values().index(MET)
        start_codon = translation_table.keys()[index_of_MET]
        index = self.sequence.find(start_codon)
        if index != -1:
            aminoacids = []
            for i in range(index, len(self.sequence) - 2, 3):
                codon = self.sequence[i:i+3]
                aminoacid = translation_table[codon]
                if aminoacid == STOP:
                    break
                
                aminoacids.append(aminoacid)
                
            return Protein(self, aminoacids)
            
    def __str__ (self):
        """ Builds a string representation of this MessengerRNA """
        return "mRNA(%s)" % str(self.sequence)


class NonCodingRNA (GRNElement):
    """ Non-coding RNA (ncRNA) """

    def __init__ (self, parent, sequence):
        """ Creates new ncRNA """
	GRNElement.__init__(self, parent)
        self.sequence = sequence

    def create_miRNAs (self, binding_size):
        """ Creates miRNAs from this ncRNA """
        return [MicroRNA(self, self.sequence[index:index + binding_size]) \
                for index in self.__find_stems(binding_size)]

    def __find_stems (self, binding_size):
        """ Checks if this ncRNA contains a stem loop """
        
        stems = []
        i = 0
        while i < len(self.sequence) - binding_size:
            segment = complement_string(self.sequence[i:i + binding_size])[::-1]
            
            # Check if the rest of the ncRNA sequence contains the reverse
            # of the complement
            index = self.sequence.find(segment, i + binding_size)
            if index != -1:
                stems.append(i)
                
                # Proceed after this stem
                # TODO: What if there is another stem in the middle of this one?
                i = index + binding_size - 1
                
            i += 1
                
        return stems
            
    def __str__ (self):
        """ Builds a string representation of this ncRNA """
        return "ncRNA(%s)" % str(self.sequence)


class MicroRNA (GRNElement):
    """ MicroRNA (miRNA) """

    def __init__ (self, parent, sequence):
        """ Creates new miRNA """
	GRNElement.__init__(self, parent)
        self.sequence = sequence

    def binds_to_mRNA(self, mRNA):
        """
        Checks if this miRNA binds to the given mRNA
        """
        return mRNA.sequence.find(complement_string(self.sequence)) != -1
        
    def __str__ (self):
        """
        Builds a string representation of this miRNA
        """
        return "miRNA(%s)" % str(self.sequence)
