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

#ifndef GENE_H_
#define GENE_H_

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include "Specie.h"

using namespace std;

class Gene : public Specie {
public:
        friend ostream& operator<<(ostream& oss, const Gene& g) {
                unsigned int i;
                
                oss << "Gene " << g.specieId << ":\n\t";
                for( i = 0; i < g.sequence.size(); i++ ) {
                        oss << g.sequence.at(i);
                        if( i != 0 && i % 100 == 0 )
                                oss << "\n\t";
                }
                oss << "(dna)\n\t";
                for( i = 0; i < g.reg_region.size(); i++ ) {
                        oss << g.reg_region.at(i);
                        if( i != 0 && i % 100 == 0 )
                                oss << "\n\t";
                }                       
                oss << "(reg region)\n\t" << g.promoter << "(promoter)\n";
                
                return oss;
        }
        
        Gene() {}
        virtual ~Gene() {}
        
        Gene( string n, string p, string d, string rg ) {
                string s( "Gene" + n );
                Specie::setID( n );
                Specie::setSequence( d );
                Specie::setParent( NULL );
                Specie::setName( s );
                
                reg_region = rg;
                active = false;
                promoter = p;
        
                //regulates = new vector<Gene*>();
        }
        
        string& getRegRegion() { return reg_region; }
        unsigned int getRegRegionSize() { return reg_region.size(); }
        bool isActive() { return active; }      
private:
        // por um vector para conter varias regioes diferentes
        string reg_region;
        string promoter;
        bool active;
};

#endif /*GENE_H_*/

