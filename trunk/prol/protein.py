from defs import *
from grn_element import *

class Protein(GRNElement):
    """ A protein. """

    def __init__ (self, parent, sequence):
	GRNElement.__init__(self, parent)
        self.sequence = sequence

    def dna_aminoacids_binding(self, func, dna, aminoacids):
        """
        Finds the binding probability between a sequence of DNA and a sequence of
        aminoacids.
        """
        return func([binding_table[(a, int(b))] / 100.0 for (b, a) in zip(dna, aminoacids)])

    def binds_to_gene (self, func, gene, binding_size, binding_threshold):
        """
        Checks if this proteins binds to some gene. 
        
        Given:
        - prot_segment is some protein segment with binding_size size
        - gene_segment is some protein dna sequence with binding_size size, from
          the gene's regulatory region
        
        The protein will bind if dna_protein_binding(...) is bigger than the given
        threshold
        
        . Parameters:
        - func: method of selection of the binding probability from all the
        probablities between bases and aminoacids (max, min, avg, ...)
        - gene: The gene.
        - binding_size: Minimum size of the binding site between the gene and the
        protein.
        - binding_threshold: Minimum value for the selected binding probability
        for binding to occur.
        """
        if len(self.sequence) >= binding_size:
            if len(gene.regulatory_region) >= binding_size:
                # for each possible protein segment with size = binding site size
                for pi in range(len(self.sequence) - binding_size + 1):
                    prot_segment = self.sequence[pi:pi + binding_size]
                    
                    # for each possible dna sequence with size = binding site size
                    for gi in range(len(gene.regulatory_region) - binding_size + 1):
                        gene_segment = gene.regulatory_region[gi:gi + binding_size]
                        
                        # Check if they bind
                        if self.dna_aminoacids_binding(func, gene_segment, prot_segment) >= binding_threshold:
                            return True
                        
        return False        

    def __str__ (self):
        """ Builds a string representation of this protein """
        return "Protein(%s)" % str(self.sequence)
