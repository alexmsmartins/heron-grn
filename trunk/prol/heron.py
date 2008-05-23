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
import gene
import rna
import protein
import graph
import random
import pydot
import pickle
import yaml
import sys
from optparse import OptionParser

# Default configuration
config = {
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

def read_configuration(filename):
    """
    Read the configuration from some file
    """
    try:
        fx = open(filename, "r")
        config = yaml.load(fx)
        fx.close()
    except IOError, (errno, strerror):
        print "Error while loading the configuration: %s" % strerror


def parse_args():
    """
    Parse command-line arguments
    """
    usage = "usage: %prog [options] <genome size> <output>"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="filename",
                      help="configuration file", metavar="FILE")
    parser.add_option("-d", "--dot", dest="dot_filename",
                      help="Write graph to the dot format", metavar="FILE")
    parser.add_option("-s", action="store_true", dest="statistics",
                      help="Print some statistics")
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
        
    if options.filename != None:
        read_configuration(options.filename)

    return (options, genome_size, output_file)


def create_graph(genes, mRNAs, ncRNAs, proteins, miRNAs, options):
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

    # Create connections between proteins and genes
    if options.verbose:
        print "   - Finding bindings between proteins and genes"
    protein_gene_bindings = 0
    for (i, protein) in zip(range(len(proteins)), proteins):
        if options.verbose:
            print "    %2.0f%%" % (100.0 * i / len(proteins)),
	    sys.stdout.flush()
        for gene in genes:
            if protein.binds_to_gene(eval(config["binding function"]), gene, \
                                     config["protein/gene binding site size"], \
                                     config["protein/gene binding threshold"]):
		wt = 1
		if random.random() <= config["protein inhibition rate"]:
		    wt = -1
                grn.add_arrow(protein, gene, wt)
		protein_gene_bindings += 1

    # Create connections between miRNAs and mRNAs
    if options.verbose:
	print
        print "   - Finding bindings between miRNAs and mRNAs"
    miRNA_mRNA_bindings = 0
    for miRNA in miRNAs:
        for mRNA in mRNAs:
            if miRNA.binds_to_mRNA(mRNA):
                grn.add_arrow(miRNA, mRNA, -1)
		miRNA_mRNA_bindings += 1

    return (grn, protein_gene_bindings, miRNA_mRNA_bindings)


if __name__ == '__main__':
    (options, genome_size, output_file) = parse_args()
    
    # Exemplo do paper
    #genes = Genome("1111121203221103300301011013232200121230320022321230302031111").get_genes()

    # Generate random genome
    if options.verbose:
        print " * Generating the random genome with size %d " % genome_size
    genes = generate_random_genome(genome_size).get_genes(config["promoter"], \
                                                          config["termination"])
    mRNAs = set()
    ncRNAs = []

    # Splice the genes
    if options.verbose:
        print " * Splicing %d genes" % len(genes)
    for gene in genes:       
        (mRNA, gene_ncRNAs) = gene.splice(config["U1 left"], config["U1 right"])
	mRNAs.add(mRNA)
	ncRNAs += gene_ncRNAs

    # Translate the mRNAs into proteins
    if options.verbose:
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

    # Create miRNAs
    if options.verbose:
        print " * Creating miRNAs from %d ncRNAs" % len(ncRNAs)
    miRNAs = []
    for ncRNA in ncRNAs:
	miRNAs += ncRNA.create_miRNAs(config["miRNA/mRNA binding site size"])

    # Create the graph
    if options.verbose:
        print " * Creating the graph"
    (grn, protein_gene_bindings, miRNA_mRNA_bindings) = \
	create_graph(genes, mRNAs, ncRNAs, proteins, miRNAs, options)

    # Save the graph to the output file
    fx = open(output_file, "w")
    pickle.dump(grn, fx)
    fx.close()

    if options.statistics:
	print
	print " * Statistics:"
        print "   - Number of genes: %d" % len(genes)
        print "   - Number of mRNAs: %d" % len(mRNAs)
        print "   - Number of proteins: %d" % len(proteins)
        print "   - Number of ncRNA: %d" % len(ncRNAs)
        print "   - Number of miRNA: %d" % len(miRNAs)
        print "   - Number of Protein/Gene bindings: %d" % protein_gene_bindings
        print "   - Number of miRNA/mRNA bindings: %d" % miRNA_mRNA_bindings

    # Save to .dot format
    if options.dot_filename != None:
        edges = []
        for node in grn.get_nodes():
            for edge in grn.nodes[node]:
                edges.append((node.id, edge.id))

        g = pydot.graph_from_edges(edges, directed=True)
        g.write(options.dot_filename);



