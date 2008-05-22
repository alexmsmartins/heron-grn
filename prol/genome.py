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


from random import randint
import re

from gene import *

def generate_random_genome (size):
    """ Generates a random genome with the given size """
    return Genome("".join([str(randint(0, 3)) for i in range(size)]))


class Genome(object):
    """ Artificial genome """

    def __init__ (self, sequence):
        """ Creates a new genome with the given sequence """
        self.sequence = sequence

    def __str__ (self):
        """ Returns a string representation of the genome """
        return str(self.sequence)

    def get_genes (self):
        """ Returns the genes written in the genome """

        promoter = "0101"
        termination = "1111"
        genes = []

        termination_re = re.compile(termination)
        promoter_re = re.compile(promoter)
        for gene_seq in termination_re.split(self.sequence):
            regions = promoter_re.split(gene_seq)
            
            # if the gene has one promoter
            if len(regions) == 2:
                regulatory, coding = regions
                genes.append(Gene(regulatory, coding))
                
        return genes
