import re
from defs import *
from protein import *

class RNA (object):
    """ Messenger RNA """

    def __init__ (self, sequence):
        """ Creates new messenger RNA """
        self.sequence = sequence
        

    def translate (self):
        """ Translate the mRNA into proteins """

        start_codon = translation_table.keys()[translation_table.values().index(MET)]
        regex = re.compile(start_codon)

        match = regex.search(self.sequence)
        if match != None:
            protein = []
            for i in range(match.start(), len(self.sequence) - 2, 3):
                codon = self.sequence[i:i+3]
                aminoacid = translation_table[codon]
                if aminoacid == STOP:
                     break

                protein.append(translation_table[codon])

        return Protein(protein)

    def __str__ (self):
        """ Builds a string representation of this RNA """
        return str(self.sequence)
