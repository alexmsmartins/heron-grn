from genome import *
from gene import *
from rna import *

genes = generate_random_genome(10000).get_genes()
(mRNA, ncRNAs) = genes[0].splice()

print mRNA.translate()


genome = Genome("1111121203221103300301011013232200121230320022321230302031111")
gene = genome.get_genes()[0]

print "Regulation region: ", gene.regulatory_region
print "Coding region: ", gene.coding_region
print "After transcription: ", gene.transcribe()

(mRNA, ncRNAs) = gene.splice()

print "ncRNAs: ", ncRNAs

print mRNA.sequence

print mRNA.translate()
