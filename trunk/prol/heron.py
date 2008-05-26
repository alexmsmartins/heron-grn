#! /usr/bin/python

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

from genome import *
from rna import NonCodingRNA
import graph
import random
import pydot
import pickle
import yaml
import sys
from optparse import OptionParser


# Default configuration
default_config = {
    "miRNA/mRNA binding site size" : 6,
    "protein/gene binding site size" : 6,
    "protein/gene binding threshold" : 0.3,
    "U1 left" : "30",
    "U1 right" : "13",
    "promoter" : "0101",
    "termination" : "1111",
    "binding function" : "lambda list: sum(list) / len(list)",
    "protein inhibition rate" : 0.5
}


class Heron:
    """
    The HeRoN model for genetic regulatory networks
    """

    def __init__(self, genome_size, config=default_config):
        """
        Initializes this instance.

        Arguments:
        - config: The model configuration with all the parameters necessary
                  to produce the graph
        """
        
        self.config = config
        self.genome_size = genome_size

    def create_network(self, window=None, event=None, verbose=False):
        """
        Creates a regulatory network from the given genome size
        """

        def advance(**args):
            """
            Event notification stuff
            """
            if window != None and event != None:
                import wx
                evt = event(**args)
                wx.PostEvent(window, evt)

        if verbose:
            print " * Generating the random genome with size %d " % self.genome_size
        genome = generate_random_genome(self.genome_size)
        genes = genome.get_genes(self.config["promoter"], \
                                     self.config["termination"])

        advance(step=1, genes=len(genes))

        mRNAs = set()
        ncRNAs = []

        # Splice the genes
        if verbose:
            print " * Splicing %d genes" % len(genes)
        for gene in genes:       
            (mRNA, gene_ncRNAs) = gene.splice(self.config["U1 left"], \
                                                  self.config["U1 right"])
            mRNAs.add(mRNA)
            ncRNAs += gene_ncRNAs

        advance(step=2)

        # Translate the mRNAs into proteins
        if verbose:
            print " * Translating %d mRNAs" % len(mRNAs)
        proteins = []
        to_remove = set()
        for mRNA in mRNAs:
            protein = mRNA.translate()
            if protein != None:
                proteins.append(protein)
            else:
                # if the mRNA does not contain the start codon, add it to the
                # list of non-coding RNA
                ncRNAs.append(NonCodingRNA(mRNA.parent, mRNA.sequence))
                to_remove.add(mRNA)

        # Remove mRNAs without the start codon from the list of mRNAs
        mRNAs -= to_remove

        advance(step=3, mRNAs=len(mRNAs), ncRNAs=len(ncRNAs), proteins=len(proteins))

        # Create miRNAs
        if verbose:
            print " * Creating miRNAs from %d ncRNAs" % len(ncRNAs)
        miRNAs = []
        for ncRNA in ncRNAs:
            miRNAs += ncRNA.create_miRNAs(self.config["miRNA/mRNA binding site size"])

        if verbose:
            print " * Found %d miRNAs" % len(miRNAs)

        advance(step=4, miRNAs=len(miRNAs))

        # Create the graph
        if verbose:
            print " * Creating the graph"
        self.grn = self.__create_graph(genes, mRNAs, ncRNAs, proteins, miRNAs, verbose, advance)

    def dump(self, output_file):
        """
        Save the graph to the output file
        """
        
        fx = open(output_file, "w")
        pickle.dump(self.grn, fx)
        fx.close()

    def save_as_dot(self, output_file):
        """
        Save the graph as a .dot file
        """
        edges = []
        for node in self.grn.get_nodes():
            for edge in self.grn.nodes[node]:
                edges.append((node.id, edge.id))
                
        g = pydot.graph_from_edges(edges, directed=True)
        g.write(output_file)

    def __create_graph(self, genes, mRNAs, ncRNAs, proteins, miRNAs, verbose, advance=lambda n:None):
        """
        Creates the graph containing the gene regulatory network
        """
        
        grn = graph.graph()
        grn.add_nodes(genes)
        grn.add_nodes(mRNAs)
        grn.add_nodes(ncRNAs)
        grn.add_nodes(proteins)
        grn.add_nodes(miRNAs)
        
        # Create connections between genes and mRNA
        for mRNA in mRNAs:
            assert isinstance(mRNA.parent, Gene)
            grn.add_arrow(mRNA.parent, mRNA, 1)
            
        # Create connections between genes and ncRNA
        for ncRNA in ncRNAs:
            assert isinstance(ncRNA.parent, Gene)
            grn.add_arrow(ncRNA.parent, ncRNA, 1)

        # Create connections between ncRNA and miRNA
        for miRNA in miRNAs:
            assert isinstance(miRNA.parent, NonCodingRNA)
            grn.add_arrow(miRNA.parent, miRNA, 1)

        # Create connections between mRNA and proteins
        for protein in proteins:
            assert isinstance(protein.parent, MessengerRNA)
            grn.add_arrow(protein.parent, protein, 1)

        # Binding function definitions, so that eval() recognizes them
        average = lambda list: sum(list) / len(list)

        # Create connections between proteins and genes
        bindings = 0
        if verbose:
            print "   - Finding bindings between proteins and genes"
        for (i, protein) in enumerate(proteins):
#            if verbose:
#                print "    %2.0f%%" % (100.0 * i / len(proteins)),
#                sys.stdout.flush()
            advance(step=4, done=i*1.0/len(proteins), protein_gene_bindings=bindings)
            for gene in genes:
                if protein.binds_to_gene(eval(self.config["binding function"]), gene, \
                                             self.config["protein/gene binding site size"], \
                                             self.config["protein/gene binding threshold"]):
                    wt = 1
                    if random.random() <= self.config["protein inhibition rate"]:
                        wt = -1
                    grn.add_arrow(protein, gene, wt)
                    bindings += 1
                    
        if verbose:
            print "     ~ Found %d protein/gene bindings" % bindings
        advance(step=5, protein_gene_bindings=bindings)
                    
        # Create connections between miRNAs and mRNAs
        if verbose:
            print "   - Finding bindings between miRNAs and mRNAs"
        bindings = 0
        for (i, miRNA) in enumerate(miRNAs):
            advance(step=5, done=i*1.0/len(miRNAs), miRNA_mRNA_bindings=bindings)
            for mRNA in mRNAs:
                if miRNA.binds_to_mRNA(mRNA):
                    grn.add_arrow(miRNA, mRNA, -1)
                    bindings += 1

        if verbose:
            print "     ~ Found %d miRNA/mRNA bindings" % bindings
        advance(step=6, done=1, miRNA_mRNA_bindings=bindings)

        return grn

       
def read_configuration(filename):
    """
    Read the configuration from some file
    """
    try:
        fx = open(filename, "r")
        config = yaml.load(fx)
        return config
    except IOError, (errno, strerror):
        print "Error while loading the configuration: %s" % strerror
    finally:
        fx.close()

def parse_args():
    """
    Parse command-line arguments
    """
    usage = "usage: %prog [options] <genome size> <output>"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="config_filename",
                      help="configuration file", metavar="FILE")
    parser.add_option("-d", "--dot", dest="dot_filename",
                      help="Write graph to the dot format", metavar="FILE")
    parser.add_option("-v", action="store_true", dest="verbose",
                      help="Verbose mode")

    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrent number of arguments")
        exit()

    try:
        genome_size = int(args[0])
    except:
        parser.error("Invalid genome size")
        exit()

    output_file = args[1]

    return (options, genome_size, output_file)

 
if __name__ == '__main__':
    (options, genome_size, output_file) = parse_args()

    config = default_config
    if options.config_filename != None:
        config = read_configuration(options.config_filename)
    
    heron = Heron(genome_size, config)
    heron.create_network(verbose=options.verbose)
    heron.dump(output_file)

    if options.dot_filename != None:
        heron.save_as_dot(options.dot_filename)    

