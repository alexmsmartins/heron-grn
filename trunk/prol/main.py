from genome import *
from gene import *
from rna import *
from protein import *
from util import *

def dna_protein_binding(func, dna, protein):
    """
    Finds the binding probability between a sequence of DNA and a sequence of
    aminoacids.
    """
    return func([binding_table[(a, int(b))] / 100.0 for (b, a) in zip(dna, protein)])


def gene_protein_binding (func, gene, protein, binding_site_size,
                          binding_threshold):
    """ Finds the binding probability between a gene and a protein. """

    if len(protein.sequence) >= binding_site_size:
        if len(gene.regulatory_region) >= binding_site_size:
            # for each possible protein segment with size = binding site size
            for pi in range(len(protein.sequence) - binding_site_size + 1):
                prot_segment = protein.sequence[pi:pi + binding_site_size]

                # for each possible dna sequence with size = binding site size
                for gi in range(len(gene.regulatory_region) - binding_site_size + 1):
                    gene_segment = gene.regulatory_region[gi:gi + binding_site_size]

                    # Check if they bind
                    if dna_protein_binding(func, gene_segment, prot_segment) >= binding_threshold:
                        return True

    return False

def find_ncRNA_stems (rna, miRNA_binding_site_size):
    """ Checks if the given ncRNA contains a stem loop """
    
    stems = []
    i = 0
    while i < len(rna) - miRNA_binding_site_size:
        segment = rna[i:i + miRNA_binding_site_size]

        # Check if the rest of the ncRNA sequence contains the reverse
        # of the complement 
        index = rna.find(complement_string(segment)[::-1],
                         i + miRNA_binding_site_size)
        if index != -1:
            stems.append(i)

            # Proceed after this stem
            # TODO: What if there is another stem in the middle of this one?
            i = index + miRNA_binding_site_size - 1

        i += 1

    return stems
        
    
                        
genes = generate_random_genome(10000).get_genes()
(mRNA, ncRNAs) = genes[0].splice()

print mRNA.translate()

# Exemplo do paper
genome = Genome("1111121203221103300301011013232200121230320022321230302031111")
gene = genome.get_genes()[0]

print "Regulation region: ", gene.regulatory_region
print "Coding region: ", gene.coding_region
print "After transcription: ", gene.transcribe()

(mRNA, ncRNAs) = gene.splice()

print "ncRNAs: ", ncRNAs

print mRNA.sequence

protein = mRNA.translate()

print protein

print gene_protein_binding(max, gene, protein, 4, 0.1)
print dna_protein_binding(lambda l: sum(l) / len(l), [T, A, G, A], [MET, VAL, VAL, ALA])

print find_ncRNA_stems("30321231133230321013000000111111222222333333", 6)
