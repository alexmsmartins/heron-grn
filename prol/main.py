from genome import *
import gene
import rna
import protein
import graph
import random

# Exemplo do paper
#genome = Genome("1111121203221103300301011013232200121230320022321230302031111")

if __name__ == '__main__':
    miRNA_binding_site_size = 6
    protein_gene_binding_site_size = 6
    protein_gene_binding_threshold = 0.3

    genes = generate_random_genome(100000).get_genes()
    mRNAs = set()
    ncRNAs = []
    for gene in genes:       
        (mRNA, gene_ncRNAs) = gene.splice()
	mRNAs.add(mRNA)
	ncRNAs += gene_ncRNAs

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

    mRNAs -= to_remove

    miRNAs = []
    for ncRNA in ncRNAs:
	miRNAs += ncRNA.create_miRNAs(miRNA_binding_site_size)

    print "Number of genes: %d" % len(genes)
    print "Number of mRNAs: %d" % len(mRNAs)
    print "Number of proteins: %d" % len(proteins)
    print "Number of ncRNA: %d" % len(ncRNAs)
    print "Number of miRNA: %d" % len(miRNAs)

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
            if protein.binds_to_gene(lambda l: sum(l) / len(l), gene, protein_gene_binding_site_size, protein_gene_binding_threshold):
                grn.add_edge(protein, gene, wt=random.choice([-1, 1]))

    # Todo
    # Create connections between miRNA and mRNA

    print len(grn.weights)
#                blah[protein] = gene

#	        print "%s binds to %s" % (str(protein), str(gene))



    # TODO
    # Find interactions between miRNAs and genes
    # Create the graph
