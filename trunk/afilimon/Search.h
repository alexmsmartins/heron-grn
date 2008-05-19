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

#include <string>
#include <iostream>
#include <assert.h>

using namespace std;

class Search {
    static const int MAXCHAR = 256;
    int d[MAXCHAR];
    int m;
    char* patt;
    int limit;
public:
    Search() { limit = 0; }
    
    int find(char *text, char *p) {
        assert(p);
            patt = p;
            m = strlen(patt);
        
            int k = 0;
        
            for( k = 0; k < MAXCHAR; k++ )
                    d[k] = m;
        
            for( k = 0; k < m - 1; k++ )
                    d[patt[k]] = m - k - 1;
                
            assert(text);
            int missed, j, i;
            int n = strlen(text);
            if( m > n )
                    return -1;
        
            k = m - 1;
        
            while( k < n ) {
                j = m - 1;
                i = k;
                missed = 0;
                //std::cout << "Starting at " << i << " in the text\n";
                while( j >= 0 && (text[i] == patt[j] || (text[i] != patt[j] && missed++ < limit)) ) {
                        //std::cout << "\tmissed = " << missed << "\n";
                        j--;
                    i--;
                }
                if( j == -1 )
                        return i + 1;
                k += d[text[k]];
            }
        
            return string::npos;
    }
    
    int find(string &text, string &p) {
        m = p.size();
            int k = 0;
        
                //cout << "Searching " << p << " in: " << text << "\n";
            for( k = 0; k < MAXCHAR; k++ )
                    d[k] = m;
        
            for( k = 0; k < m - 1; k++ )
                    d[p[k]] = m - k - 1;
                
            int missed, j, i;
            int n = text.size();
            if( m > n )
                    return -1;
        
            k = m - 1;
                
            while( k < n ) {
                j = m - 1;
                i = k;
                missed = 0;
                //std::cout << "Starting at " << i << " in the text\n";
                while( j >= 0 && (text[i] == p[j] || (text[i] != p[j] && missed++ < limit)) ) {
                        //std::cout << "\tmissed = " << missed << "\n";
                        j--;
                    i--;
                }
                if( j == -1 )
                        return i + 1;
                k += d[text[k]];
            }
        
            return string::npos;
        }
                
        void setLimit(int l) { limit = l; }
};

