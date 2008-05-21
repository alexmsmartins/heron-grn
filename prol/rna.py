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
        start_codon = translation_table.keys()[translation_table.values().index(MET)]
        index = self.sequence.find(start_codon)
        if index != -1:
            aminoacids = []
            for i in range(index, len(self.sequence) - 2, 3):
                codon = self.sequence[i:i+3]
                aminoacid = translation_table[codon]
                if aminoacid == STOP:
                    break
                
                aminoacids.append(translation_table[codon])
                
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
        return [MicroRNA(self, self.sequence[index:index + binding_size]) for index in self.__find_stems(binding_size)]

    def __find_stems (self, binding_size):
        """ Checks if this ncRNA contains a stem loop """
        
        stems = []
        i = 0
        while i < len(self.sequence) - binding_size:
            segment = self.sequence[i:i + binding_size]
            
            # Check if the rest of the ncRNA sequence contains the reverse
            # of the complement 
            index = self.sequence.find(complement_string(segment)[::-1], i + binding_size)
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
        
    def __str__ (self):
        """ Builds a string representation of this miRNA """
        return "miRNA(%s)" % str(self.sequence)
