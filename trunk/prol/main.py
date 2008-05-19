from genome import *
from gene import *
from rna import *
from protein import *

def dna_protein_binding(func, dna, protein):
    """
    Finds the binding probability between a sequence of DNA and a sequence of
    aminoacids.
    """
    return func([binding_table[(a, b)] / 100.0 for (b, a) in zip(dna, protein)])

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

print mRNA.translate()

print dna_protein_binding(lambda l: sum(l) / len(l), [T, A, G, A], [MET, VAL, VAL, ALA])
