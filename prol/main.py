from genome import *
import gene
import rna
import protein

# Exemplo do paper
#genome = Genome("1111121203221103300301011013232200121230320022321230302031111")

if __name__ == '__main__':
    miRNA_binding_site_size = 6
    protein_gene_binding_site_size = 6
    protein_gene_binding_threshold = 0.3

    genes = generate_random_genome(100000).get_genes()
    mRNAs = []
    ncRNAs = []
    for gene in genes:
        (mRNA, gene_ncRNAs) = gene.splice()
	mRNAs.append(mRNA)
	ncRNAs += gene_ncRNAs

    proteins = []
    for mRNA in mRNAs:
        protein = mRNA.translate()
	if protein != None:
	    proteins.append(protein)
	else:
	    # if the mRNA does not contain the start codon, add it to the
	    # list of non-coding RNA
	    ncRNAs.append(NonCodingRNA(mRNA.parent, mRNA.sequence))

    miRNAs = []
    for ncRNA in ncRNAs:
	miRNAs += ncRNA.create_miRNAs(miRNA_binding_site_size)

    print "Number of genes: %d" % len(genes)
    print "Number of proteins: %d" % len(proteins)
    print "Number of ncRNA: %d" % len(ncRNAs)
    print "Number of miRNA: %d" % len(miRNAs)
        
    for protein in proteins:
	for gene in genes:
	    if protein.binds_to_gene(max, gene, protein_gene_binding_site_size, protein_gene_binding_threshold):
	        print "%s binds to %s" % (str(protein), str(gene))

    # TODO
    # Find interactions between miRNAs and genes
    # Create the graph
