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


# Bases
T, A, G, C = range(4)

# Uracil replaces thymine in RNA
U = T

# Aminoacids
ARG, LYS, HIS, GLU, GLN, ASP, ASN, TYR, PHE, TRP, LEU, ILE, VAL, MET, THR, \
SER, CYS, ALA, GLY, PRO, STOP = range(21)

# This table relates codons with aminoacids for translation
# This is only an auxiliary table where the keys are tuples. The table that is
# really used is translation_table and is defined after this table.
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

# Similar to the previous table, but the element key is a string of bases
# instead of a tuple (useful because the RNA sequence is also a string)
translation_table = {}
for key in aux_table.keys():
    translation_table["".join(map(str, key))] = aux_table[key]

# Binding probabilities between aminoacids and bases
binding_table = { 
    # Arginine (R)
    (ARG, A) : 19.6,
    (ARG, T) : 12.2,
    (ARG, C) : 24.1,
    (ARG, G) : 35.7,

    # Lysine (K)
    (LYS, A) : 23.7,
    (LYS, T) : 16.3,
    (LYS, C) : 22.8,
    (LYS, G) : 30.7,

    # Histidine (H)
    (HIS, A) : 25.3,
    (HIS, T) : 14.2,
    (HIS, C) : 16.2,
    (HIS, G) : 37.7,

    # Glutamic acid (E)
    (GLU, A) : 19.1,
    (GLU, T) : 4.8,
    (GLU, C) : 34.8,
    (GLU, G) : 33.0,

    # Glutamine (Q)
    (GLN, A) : 28,
    (GLN, T) : 13.7,
    (GLN, C) : 17.7,
    (GLN, G) : 29.4,

    # Aspartic acid (D)
    (ASP, A) : 13.3,
    (ASP, T) : 1.5,
    (ASP, C) : 34.2,
    (ASP, G) : 37,

    # Asparagine (N)
    (ASN, A) : 25.5,
    (ASN, T) : 17.7,
    (ASN, C) : 20,
    (ASN, G) : 23.9,
    
    # Tyrosine (Y) 
    (TYR, A) : 28.4,
    (TYR, T) : 15,
    (TYR, C) : 27.4,
    (TYR, G) : 23.6,
    
    # Phenylalanine (F)
    (PHE, A) : 17.7,
    (PHE, T) : 17.7,
    (PHE, C) : 24.1,
    (PHE, G) : 40.5,
    
    # Trytophan (W)
    (TRP, A) : 14.4,
    (TRP, T) : 21.8,
    (TRP, C) : 30.2,
    (TRP, G) : 24.8,
    
    # Leucine (L) 
    (LEU, A) : 9.5,
    (LEU, T) : 19.4,
    (LEU, C) : 31.1,
    (LEU, G) : 30.2,
    
    # Isoleucine (I)
    (ILE, A) : 21.4,
    (ILE, T) : 11.4,
    (ILE, C) : 26.4,
    (ILE, G) : 30.8,
    
    # Valine (V)
    (VAL, A) : 25,
    (VAL, T) : 13.3,
    (VAL, C) : 35.3,
    (VAL, G) : 20,
    
    # Methionine (M)
    (MET, A) : 22.1,
    (MET, T) : 9.8,
    (MET, C) : 27.9,
    (MET, G) : 22.1,
    
    # Threonine (T)
    (THR, A) : 24.6,
    (THR, T) : 23.1,
    (THR, C) : 20.2,
    (THR, G) : 27.8,
    
    # Serine (S)
    (SER, A) : 28.2,
    (SER, T) : 19.7,
    (SER, C) : 20.9,
    (SER, G) : 27.2,
    
    # Cysteine (C)
    (CYS, A) : 29.1,
    (CYS, T) : 23.1,
    (CYS, C) : 18.8,
    (CYS, G) : 24.8,
    
    # Alanine (A)
    (ALA, A) : 24.2,
    (ALA, T) : 24.6,
    (ALA, C) : 17.3,
    (ALA, G) : 24,
    
    # Glycine (G)
    (GLY, A) : 20.1,
    (GLY, T) : 17,
    (GLY, C) : 22.9,
    (GLY, G) : 32.1,
    
    # Proline (P)
    (PRO, A) : 37,
    (PRO, T) : 2,
    (PRO, C) : 11,
    (PRO, G) : 21
}
