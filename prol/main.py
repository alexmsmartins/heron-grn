#! /usr/bin/python
# -*- coding: utf-8 -*-

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

avg = lambda l: sum(l) / len(l)
miRNA_binding_site_size = 6
protein_gene_binding_site_size = 6
protein_gene_binding_threshold = 0.34


def create_graph(genes, mRNAs, ncRNAs, proteins, miRNAs):
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
        grn.add_edge(mRNA.parent, mRNA, 1)

    # Create connections between genes and ncRNA
    for ncRNA in ncRNAs:
        assert isinstance(ncRNA.parent, Gene)
        grn.add_edge(ncRNA.parent, ncRNA, 1)

    # Create connections between ncRNA and miRNA
    for miRNA in miRNAs:
        assert isinstance(miRNA.parent, NonCodingRNA)
        grn.add_edge(miRNA.parent, miRNA, 1)

    # Create connections between mRNA and proteins
    for protein in proteins:
        assert isinstance(protein.parent, MessengerRNA)
        grn.add_edge(protein.parent, protein, 1)

    # Create connections between proteins and genes
    for protein in proteins:
        for gene in genes:
            if protein.binds_to_gene(avg, gene, \
                                     protein_gene_binding_site_size, \
                                     protein_gene_binding_threshold):
                grn.add_edge(protein, gene, wt=random.choice([-1, 1]))

    # Create connections between miRNAs and mRNAs
    for miRNA in miRNAs:
        for mRNA in mRNAs:
            if miRNA.binds_to_mRNA(mRNA):
                grn.add_edge(miRNA, mRNA, -1)

    return grn


if __name__ == '__main__':
    # Exemplo do paper
    #genes = Genome("1111121203221103300301011013232200121230320022321230302031111").get_genes()

    # Generate random genome
    genes = generate_random_genome(100000).get_genes()
    mRNAs = set()
    ncRNAs = []

    # Splice the genes
    for gene in genes:       
        (mRNA, gene_ncRNAs) = gene.splice()
	mRNAs.add(mRNA)
	ncRNAs += gene_ncRNAs

    # Translate the mRNAs into proteins
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
    miRNAs = []
    for ncRNA in ncRNAs:
	miRNAs += ncRNA.create_miRNAs(miRNA_binding_site_size)

    # Create the graph
    grn = create_graph(genes, mRNAs, ncRNAs, proteins, miRNAs)

    edges = []
    for i in range(len(grn.get_nodes())):
        node = grn.get_nodes()[i]
        for edge in grn.nodes[node]:
            edges += (str(node), str(edge))

#    for node in grn.get_nodes():
#        for edge in grn.nodes[node]:
#            edges += (str(node), str(edge))
    
    g = pydot.graph_from_edges(edges)
    g.write('banana.dot');


    print "Number of genes: %d" % len(genes)
    print "Number of mRNAs: %d" % len(mRNAs)
    print "Number of proteins: %d" % len(proteins)
    print "Number of ncRNA: %d" % len(ncRNAs)
    print "Number of miRNA: %d" % len(miRNAs)
