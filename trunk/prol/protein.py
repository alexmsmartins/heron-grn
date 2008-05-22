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


from defs import *
from grn_element import *

class Protein(GRNElement):
    """ A protein. """

    def __init__ (self, parent, sequence):
	GRNElement.__init__(self, parent)
        self.sequence = sequence

    def dna_aminoacids_bindings(self, dna, aminoacids):
        """
        Returns the binding probabilities between a sequence of DNA
        and a sequence of aminoacids.
        """
        return [binding_table[(a, int(b))] / 100.0
                for (b, a) in zip(dna, aminoacids)]

    def binds_to_gene (self, func, gene, binding_size, binding_threshold):
        """
        Checks if this proteins binds to some gene. 
        
        Given:
        - prot_seg is some protein segment with binding_size size
        - gene_seg is some protein dna sequence with binding_size size, from
          the gene's regulatory region
        
        The protein will bind if dna_protein_binding(...) is bigger than the
        given threshold
        
        . Parameters:
        - func: method of selection of the binding probability from all the
          probablities between bases and aminoacids (max, min, avg, ...)
        - gene: The gene.
        - binding_size: Minimum size of the binding site between the gene and
          the protein.
        - binding_threshold: Minimum value for the selected binding probability
          for binding to occur.
        """
        if len(self.sequence) >= binding_size:
            if len(gene.regulatory_region) >= binding_size:
                # foreach possible protein segment with size = binding site size
                for pi in range(len(self.sequence) - binding_size + 1):
                    prot_seg = self.sequence[pi:pi + binding_size]
                    
                    # foreach possible dna seq. with size = binding site size
                    for gi in range(len(gene.regulatory_region) - \
                                    binding_size + 1):
                        gene_seg = gene.regulatory_region[gi:gi + binding_size]
                        
                        # Check if they bind
                        if func(self.dna_aminoacids_bindings(gene_seg, prot_seg)) \
                               >= binding_threshold:
                            return True
                        
        return False        

    def __str__ (self):
        """ Builds a string representation of this protein """
        return "Protein(%s)" % str(self.sequence)
