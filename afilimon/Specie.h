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

#ifndef SPECIE_H_
#define SPECIE_H_

#include <string>

using namespace std;

class Specie {
public:
        Specie() {
                concentration = 0;
                parent = NULL;
        }
        
        Specie( string n, string s ) {
                setID( n );
                setSequence( "" );
                setParent( NULL );
                setName( s );
        }
        
        virtual ~Specie() { delete parent; }

        friend ostream& operator<<(ostream& oss, const Specie& p) {
                oss << "Specie " << p.name << ":\n\t";
                for( unsigned int i = 0; i < p.sequence.size(); i++ ) {
                        oss << p.sequence.at(i);
                        if( i != 0 && i % 100 == 0 )
                                oss << "\n\t";
                }
                oss << "(sequence)\n\t";
                if( p.parent != NULL )
                        oss << p.parent->getName() << "(parent specie) ";
                oss << "\n";
                return oss;
        }
        
        friend bool operator<(const Specie& s, const Specie& p) {
                return s.specieId == p.specieId;
        }
        
        void setID( string i ) { specieId = i; }
        void setName( string n ) { name = n; }
        void setSequence( string s ) { sequence = s; }
        void setParent( Specie* p ) { parent = p; }
        void setConcentration( float n ) { concentration = n; }
        
        string getID() { return specieId; }
        string getName() { return name; }
        string& getSequence() { return sequence; }
        float getConcentration() { return concentration; }
        Specie* getParent() { return parent; }
protected:
        string specieId;
        string name;
        string sequence;
        float concentration;
        
        /* The specie object that produced this one. 
         * NULL if there is no parent. */
        Specie* parent;
};

#endif /*SPECIE_H_*/

