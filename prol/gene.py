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


import re
from grn_element import *
from rna import *
from util import *

class Gene(GRNElement):
    """ Just a simple Gene """

    def __init__ (self, regulatory_region, coding_region):
        """ 
        Initializes a new gene with the given regulatory and coding 
        regions 
        """
        self.regulatory_region = regulatory_region
        self.coding_region = coding_region

    def transcribe (self):
        """ Transcribes the coding region into precursor mRNA """

        return "".join([str(complement(int(x))) for x in self.coding_region])

    def splice (self):
        """ Splices the pre-mRNA into mRNA and introns """
        
        U1_left = "30"
        U1_right = "13"

        intron_regex = re.compile(U1_left + ".+?" + U1_right)
        pre_mRNA = self.transcribe()

        introns = [NonCodingRNA(self, match) \
                   for match in intron_regex.findall(pre_mRNA)]
        mRNA = "".join(intron_regex.split(pre_mRNA))

        return (MessengerRNA(self, mRNA), introns)       

    def __str__ (self):
        """ Builds a string representation of this gene """
        return "Gene(%s, %s)" % (str(self.regulatory_region),
                                 str(self.coding_region))
        
        
    
