import re

T, A, G, C = range(4)
U = T
ARG, LYS, HIS, GLU, GLN, ASP, ASN, TYR, PHE, TRP, LEU, ILE, VAL, MET, THR, SER, CYS, ALA, GLY, PRO, STOP = range(21)

aux_table = {
    (T, T, T) : PHE,
    (T, T, C) : PHE,
    (T, T, A) : LEU,
    (T, T, G) : LEU,

    (T, A, T) : TYR,
    (T, A, C) : TYR,
    (T, A, A) : STOP,
    (T, A, G) : STOP,

    (T, G, T) : CYS,
    (T, G, C) : CYS,
    (T, G, A) : STOP,
    (T, G, G) : TRP,

    (C, A, T) : HIS,
    (C, A, C) : HIS,
    (C, A, A) : GLN,
    (C, A, G) : GLN,

    (A, T, T) : ILE,
    (A, T, C) : ILE,
    (A, T, A) : ILE,
    (A, T, G) : MET,

    (A, A, T) : ASN,
    (A, A, C) : ASN,
    (A, A, A) : LYS,
    (A, A, G) : LYS,

    (A, G, T) : SER,
    (A, G, C) : SER,
    (A, G, A) : ARG,
    (A, G, G) : ARG,

    (G, A, T) : ASP,
    (G, A, C) : ASP,
    (G, A ,A) : GLU,
    (G, A ,G) : GLU,
}

for i in range(4):
    aux_table[(T, C, i)] = SER
    aux_table[(C, T, i)] = LEU
    aux_table[(C, C, i)] = PRO
    aux_table[(C, G, i)] = ARG
    aux_table[(A, C, i)] = THR
    aux_table[(G, T, i)] = VAL
    aux_table[(G, C, i)] = ALA
    aux_table[(G, G, i)] = GLY

translation_table = {}
for key in aux_table.keys():
    translation_table["".join(map(str, key))] = aux_table[key]

class RNA (object):
    """ Messenger RNA """

    def __init__ (self, sequence):
        """ Creates new messenger RNA """
        self.sequence = sequence
        

    def translate (self):
        """ Translate the mRNA into proteins """

        start_codon = translation_table.keys()[translation_table.values().index(MET)]
        regex = re.compile(start_codon)

        match = regex.search(self.sequence)
        if match != None:
            protein = []
            for i in range(match.start(), len(self.sequence) - 2, 3):
                codon = self.sequence[i:i+3]
                aminoacid = translation_table[codon]
                if aminoacid == STOP:
                     break

                protein.append(translation_table[codon])

        return protein
