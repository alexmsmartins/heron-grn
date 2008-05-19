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

#ifndef IDGENERATOR_H_
#define IDGENERATOR_H_

class IdGenerator {
public:
        IdGenerator() {
                id = 0;
        }
        
        virtual ~IdGenerator() {}

        string getCurrentID() {
                char v[7];
                sprintf( v, "%d", id );
                id++;
                
                return string(v); 
        }
        
        void startAt( int i ) {
                id = i;
        }
private:
        int id;
};

#endif /*IDGENERATOR_H_*/

