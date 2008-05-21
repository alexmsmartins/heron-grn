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

        introns = [NonCodingRNA(self, match) for match in intron_regex.findall(pre_mRNA)]
        mRNA = "".join(intron_regex.split(pre_mRNA))

        return (MessengerRNA(self, mRNA), introns)       

    def __str__ (self):
        """ Builds a string representation of this gene """
        return "Gene(%s, %s)" % (str(self.regulatory_region),
                                 str(self.coding_region))
        
        
    
