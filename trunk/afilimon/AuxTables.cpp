/*
        Copyright 2007 Angela Goncalves
        
        This file is part of HeRoN.

    HeRoN is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    HeRoN is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

int translation_table[4][4][4];
float binding_table[20][4];
char one_letter_protein_encoding[21];

// bases
enum { T, A, G, C };
enum { U };
// proteinas
enum { ARG, LYS, HIS, GLU, GLN, ASP, ASN, TYR, PHE, TRP, LEU, ILE, VAL, MET, THR, SER, CYS, ALA, GLY, PRO, STOP };

void initialize_translation_table() {
        
        translation_table[T][T][T] = PHE;
        translation_table[T][T][C] = PHE;
        translation_table[T][T][A] = LEU;
        translation_table[T][T][G] = LEU;
        
        translation_table[T][A][T] = TYR;
        translation_table[T][A][C] = TYR;
        translation_table[T][A][A] = STOP;
        translation_table[T][A][G] = STOP;
        
        translation_table[T][G][T] = CYS;
        translation_table[T][G][C] = CYS;
        translation_table[T][G][A] = STOP;
        translation_table[T][G][G] = TRP;
        
        translation_table[C][A][T] = HIS;
        translation_table[C][A][C] = HIS;
        translation_table[C][A][A] = GLN;
        translation_table[C][A][G] = GLN;
        
        translation_table[A][T][T] = ILE;
        translation_table[A][T][C] = ILE;
        translation_table[A][T][A] = ILE;
        translation_table[A][T][G] = MET;
        
        translation_table[A][A][T] = ASN;
        translation_table[A][A][C] = ASN;
        translation_table[A][A][A] = LYS;
        translation_table[A][A][G] = LYS;
        
        translation_table[A][G][T] = SER;
        translation_table[A][G][C] = SER;
        translation_table[A][G][A] = ARG;
        translation_table[A][G][G] = ARG;
        
        translation_table[G][A][T] = ASP;
        translation_table[G][A][C] = ASP;
        translation_table[G][A][A] = GLU;
        translation_table[G][A][G] = GLU;

        for( int i = 0; i < 4; i++ ) {
                translation_table[T][C][i] = SER;
                translation_table[C][T][i] = LEU;
                translation_table[C][C][i] = PRO;
                translation_table[C][G][i] = ARG;
                translation_table[A][C][i] = THR;
                translation_table[G][T][i] = VAL;
                translation_table[G][C][i] = ALA;
                translation_table[G][G][i] = GLY;
        }
}

void initialize_one_letter_protein_encoding_table() {
        one_letter_protein_encoding[ARG] = 'R';
        one_letter_protein_encoding[LYS] = 'K';
        one_letter_protein_encoding[HIS] = 'H';
        one_letter_protein_encoding[GLU] = 'E';
        one_letter_protein_encoding[GLN] = 'Q';
        one_letter_protein_encoding[ASP] = 'D';
        one_letter_protein_encoding[ASN] = 'N';
        one_letter_protein_encoding[TYR] = 'Y';
        one_letter_protein_encoding[PHE] = 'F';
        one_letter_protein_encoding[TRP] = 'W';
        one_letter_protein_encoding[LEU] = 'L';
        one_letter_protein_encoding[ILE] = 'I';
        one_letter_protein_encoding[VAL] = 'V';
        one_letter_protein_encoding[MET] = 'M';
        one_letter_protein_encoding[THR] = 'T';
        one_letter_protein_encoding[SER] = 'S';
        one_letter_protein_encoding[CYS] = 'C';
        one_letter_protein_encoding[ALA] = 'A';
        one_letter_protein_encoding[GLY] = 'G';
        one_letter_protein_encoding[PRO] = 'P';
        one_letter_protein_encoding[STOP] = '-';
}

void initialize_binding_table() {
        // Arginine (R)
        binding_table[ARG][A] = 19.6;
        binding_table[ARG][T] = 12.2;
        binding_table[ARG][C] = 24.1; // 0
        binding_table[ARG][G] = 35.7; // 3
        
        // Lysine (K)
        binding_table[LYS][A] = 23.7;
        binding_table[LYS][T] = 16.3; // 3
        binding_table[LYS][C] = 22.8; // 0
        binding_table[LYS][G] = 30.7; // 3
        
        // Histidine (H)
        binding_table[HIS][A] = 25.3;
        binding_table[HIS][T] = 14.2;
        binding_table[HIS][C] = 16.2;
        binding_table[HIS][G] = 37.7;
        
        // Glutamic acid (E)
        binding_table[GLU][A] = 19.1;
        binding_table[GLU][T] = 4.8;
        binding_table[GLU][C] = 34.8;
        binding_table[GLU][G] = 33.0; // 0

        // Glutamine (Q)
        binding_table[GLN][A] = 28;
        binding_table[GLN][T] = 13.7;
        binding_table[GLN][C] = 17.7;
        binding_table[GLN][G] = 29.4;

        // Aspartic acid (D)
        binding_table[ASP][A] = 13.3;
        binding_table[ASP][T] = 1.5; // 0
        binding_table[ASP][C] = 34.2;
        binding_table[ASP][G] = 37; // 0
        
        // Asparagine (N)
        binding_table[ASN][A] = 25.5;
        binding_table[ASN][T] = 17.7;
        binding_table[ASN][C] = 20;
        binding_table[ASN][G] = 23.9;
        
        // Tyrosine (Y) 
        binding_table[TYR][A] = 28.4;
        binding_table[TYR][T] = 15; // 3
        binding_table[TYR][C] = 27.4; // 2
        binding_table[TYR][G] = 23.6;
        
        // Phenylalanine (F)
        binding_table[PHE][A] = 17.7;
        binding_table[PHE][T] = 17.7; // 3
        binding_table[PHE][C] = 24.1; // 2
        binding_table[PHE][G] = 40.5;
        
        // Trytophan (W)
        binding_table[TRP][A] = 14.4;
        binding_table[TRP][T] = 21.8; // 3
        binding_table[TRP][C] = 30.2; // 2
        binding_table[TRP][G] = 24.8;
        
        // Leucine (L) 
        binding_table[LEU][A] = 9.5;
        binding_table[LEU][T] = 19.4; // 3
        binding_table[LEU][C] = 31.1; // 2
        binding_table[LEU][G] = 30.2;
        
        // Isoleucine (I)
        binding_table[ILE][A] = 21.4;
        binding_table[ILE][T] = 11.4; // 3
        binding_table[ILE][C] = 26.4; // 2
        binding_table[ILE][G] = 30.8;

        // Valine (V)
        binding_table[VAL][A] = 25;
        binding_table[VAL][T] = 13.3; // 3
        binding_table[VAL][C] = 35.3; // 2
        binding_table[VAL][G] = 20;
        
        // Methionine (M)
        binding_table[MET][A] = 22.1;
        binding_table[MET][T] = 9.8; // 3
        binding_table[MET][C] = 27.9; // 2
        binding_table[MET][G] = 22.1;
        
        // Threonine (T)
        binding_table[THR][A] = 24.6;
        binding_table[THR][T] = 23.1; // 3
        binding_table[THR][C] = 20.2; // 2
        binding_table[THR][G] = 27.8;
        
        // Serine (S)
        binding_table[SER][A] = 28.2;
        binding_table[SER][T] = 19.7;
        binding_table[SER][C] = 20.9;
        binding_table[SER][G] = 27.2;
        
        // Cysteine (C)
        binding_table[CYS][A] = 29.1;
        binding_table[CYS][T] = 23.1;
        binding_table[CYS][C] = 18.8;
        binding_table[CYS][G] = 24.8;
        
        // Alanine (A)
        binding_table[ALA][A] = 24.2;
        binding_table[ALA][T] = 24.6; // 3
        binding_table[ALA][C] = 17.3; // 0
        binding_table[ALA][G] = 24; // 0
        
        // Glycine (G)
        binding_table[GLY][A] = 20.1;
        binding_table[GLY][T] = 17;
        binding_table[GLY][C] = 22.9;
        binding_table[GLY][G] = 32.1;

        // Proline (P)
        binding_table[PRO][A] = 37;
        binding_table[PRO][T] = 2;
        binding_table[PRO][C] = 11;
        binding_table[PRO][G] = 21;
}

