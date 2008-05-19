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

#ifndef PROTEIN_H_
#define PROTEIN_H_

#include "Specie.h"

class Protein : public Specie {
public:
        friend ostream& operator<<(ostream& oss, const Protein& p) {
                oss << "Protein " << p.specieId << ": " << p.sequence << "(sequence) " 
                    << p.parent->getID() << "(parent gene) " << "\n";
                return oss;
        }
        
        Protein() {}

        Protein( string id, Specie* g, string s ) {
                string n( "Protein" + id );
                Specie::setID( id );
                Specie::setParent( g );
                Specie::setSequence( s );
                Specie::setName( n );
                dummy = 0;
        }
        
        virtual ~Protein() {}
        
        int getSignal( int i ) {
                /*if( i > 0 )
                        if( sequence[sequence.size()-1] == '0' ||
                                i > 1 && sequence[sequence.size()-1] == '1' ||
                                i > 2 && sequence[sequence.size()-1] == '2' ||
                                i > 3 && sequence[sequence.size()-1] == '3' )
                                return -1;
                */
                if( i > 0 && sequence[sequence.size()-1] < '0'+i )
                        return -1;
                return 1;
        }
        
        vector<int>* getIntSeq() { return intSeq; }
        void setIntSeq( vector<int>* v ) { intSeq = v; }
private:
        int dummy;
        vector<int>* intSeq;
};

#endif /*PROTEIN_H_*/

