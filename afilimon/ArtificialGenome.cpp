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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <map>
#include <set>

#include "Specie.h"
#include "Gene.h"
#include "Protein.h"
#include "AGraph.h"
#include "MersenneTwister.h"
//#include "wildcards.hh"
//#include "MatchSubstring.h"
#include "Search.h"
#include "IdGenerator.h"
#include "AuxTables.cpp"

using namespace std;
using std::setw;

vector<string> _promoter;
unsigned int _gene_size;
unsigned int _binding_site_size;
unsigned int _miRNA_binding_site_size;
short int _num_bases;
short int _binding_site;
short int _algorithm;
int _genome_size;
int _inhibition;
int _num_steps;
int _plot_point_pt;
string _genome_file;
string _graph_file;
string _size_file;
string _input_con_file;
string _output_con_file;
string _inhibition_type;
string _termination_seq;
string _stat_file;
string _binding_choice;
string _plot_file;
string _plot_point_type;
string _save_network_file;
string _read_network_file;
string _left_U1;
string _right_U1;
float _U1_match;
float _promoter_match;
float _inhibition_rate;
float _active_genes;
float _binding_threshold;
double _active_genes_normal_dist_mean;
double _active_genes_normal_dist_var;
char _active_genes_type;

MTRand mtrand;
IdGenerator idgen;

// algorithms
enum { ARTGEN, ONE_LEVEL, TWO_LEVELS };

float average( int number_list[], int size ) {
        float sum = 0;
        
        for( int i = 0; i < size; i++ )
                sum += number_list[i];
        return sum / size;
}

int hamdist( string& s1, string& s2 ) { 
        int distance = 0;  
   
        if( s1.size() != s2.size() ) {
                cout << "Error calculation hamming distance.\n";
                return -1;
        }
        
        for( unsigned int i = 0; i < s1.size(); i++ )
        if( s1[i] != s2[i] )
                distance++;
        return distance;
}

string toProteinEncoding( vector<int>* p ) {
        char c;
        string s( "" );
        
        for( unsigned int i = 0; i < p->size(); i++ ) { 
                //cout << "converting " << p->at(i) << " of a vector with size: " << p->size() << endl;
                c = one_letter_protein_encoding[p->at(i)];
                s += c;
        }
        return s;
}

int* getProteinIntValue( string& p, int pstart, int psize ) {
        int* nv = new int[psize];
        
        for( int i = 0; i < psize; i++ ) {
                switch( p[pstart+i] ) {
                        case 'R':
                                nv[i] = ARG;
                                break;
                        case 'K':
                                nv[i] = LYS;
                                break;
                        case 'H':
                                nv[i] = HIS;
                                break;
                        case 'E':
                                nv[i] = GLU;
                                break;
                        case 'Q':
                                nv[i] = GLN;
                                break;
                        case 'D':
                                nv[i] = ASP;
                                break;
                        case 'N':
                                nv[i] = ASN;
                                break;
                        case 'Y':
                                nv[i] = TYR;
                                break;
                        case 'F':
                                nv[i] = PHE;
                                break;
                        case 'W':
                                nv[i] = TRP;
                                break;
                        case 'L':
                                nv[i] = LEU;
                                break;
                        case 'I':
                                nv[i] = ILE;
                                break;
                        case 'V':
                                nv[i] = VAL;
                                break;
                        case 'M':
                                nv[i] = MET;
                                break;
                        case 'T':
                                nv[i] = THR;
                                break;
                        case 'S':
                                nv[i] = SER;
                                break;
                        case 'C':
                                nv[i] = CYS;
                                break;
                        case 'A':
                                nv[i] = ALA;
                                break;
                        case 'G':
                                nv[i] = GLY;
                                break;
                        case 'P':
                                nv[i] = PRO;
                                break;
                        case '-':
                                nv[i] = STOP;
                                break;
                        default:
                                nv[i] = -1;
                }
        }
        return nv;
}

float dnaProteinBinding( string& p, unsigned int pstart, int psize, string& d, unsigned int dstart ) {
        float value = 0;
        int* ip = getProteinIntValue( p, pstart, psize );

        if( _binding_choice == "avg" ) {
                for( int i = 0; i < psize; i++ ) 
                        value += binding_table[ip[i]][d[dstart+i]-48];
                value = value / psize;
        } else if( _binding_choice == "max" ) {
                for( int i = 0; i < psize; i++ ) {
                        if( binding_table[ip[i]][d[dstart+i]-48] > value )
                                value = binding_table[ip[i]][d[dstart+i]-48];
                }
        } else if( _binding_choice == "min" ) {
                value = binding_table[ip[0]][d[dstart]-48];
                for( int i = 1; i < psize; i++ ) {
                        if( binding_table[ip[i]][d[dstart+i]-48] < value )
                                value = binding_table[ip[i]][d[dstart+i]-48];
                }
        } else if( _binding_choice == "random" ) {
                int r = mtrand.randInt(psize-1);
                value = binding_table[ip[r]][d[dstart+r]-48];
        }
        
        delete[] ip;
        return value;
}

void writeNetworkPlotSettings() {
        ofstream setfile( (_plot_file+".input").c_str() );
        setfile << "set title \"Network\"\n";
        setfile << "set xlabel \"Steps\"\n";
        setfile << "set ylabel \"Specie\"\n";
        setfile << "set pointsize 0.2\n";
        setfile << "set key outside\n";
        setfile << "load '" << _plot_file+".input2'\n";
        setfile << "set terminal png size 1000,600\n";
        setfile << "set output '" << _plot_file+".png'\n";
        setfile << "replot\n";
        setfile.close();
}

void writeStatisticsPlotSettingsToFile( string file, string xlabel, string ylabel ) {
        ofstream setfile( (file+".input").c_str() );
        setfile << "set xlabel \"" << xlabel << "\"\n";
        setfile << "set ylabel \"" << ylabel << "\"\n";
        setfile << "set boxwidth 1\n";
        setfile << "set style fill solid border -1\n";
        setfile << "set terminal png size 600,360\n";
        setfile << "set output '" << file+".png'\n";
        setfile << "set xrange[-1:]\n";
        setfile << "set yrange[0:]\n";
        setfile << "plot '" << file << ".dat' with boxes";
        setfile.close();
}

void executeNetwork( AGraph* graph ) {
        vector<Specie*> on, off;
        AVertex* v;
        Specie* s1;
        Specie* s2;
        int n_genes = 0;
        
        for( graph->it_first(); !graph->it_isDone(); graph->it_next() ) {
                //printf( "%2s ", graph->it_current()->getData()->getID().c_str() );
                n_genes++;
        }
        
        writeNetworkPlotSettings();
        
        ofstream file( (_plot_file+".input2").c_str() ); // ios::app 
        file << "set xrange[0:" << _num_steps+1 << "]\n";
        file << "set yrange[-1:" << n_genes << "]\n";
        file << "plot '" << _plot_file << ".dat' with " << _plot_point_type;
        if( _plot_point_pt > 0 )
                file << " pt " << _plot_point_pt;
        file.close();
        
        ofstream data( (_plot_file+".dat").c_str() );
        
        for( int i = 0; i < _num_steps; i++ ) {
                n_genes = 0;
                on.clear();
                off.clear();
                
                for( graph->it_first(); !graph->it_isDone(); graph->it_next() ) {
                        v = graph->it_current();
                        s1 = v->getData();
                        
                        //cout << s1->getName() << ":\n";
                        
                        if( s1->getConcentration() == 1 ) {
                                data << i << "\t" << s1->getID() << "\n";
                                                        
                                for( v->it_first(); !v->it_isDone(); v->it_next() ) {
                                        s2 = v->it_current();
                                
                                        if( v->currentWeight() >= 0 ) {
                                                on.push_back( s2 );
                                        } else { 
                                                off.push_back( s2 );
                                        }
                                }
                                
                                /* 1. Any gene not regulated by a currently
                                 * active gene is turned off.
                                 */
                                s1->setConcentration( 0 );
                        }
                                
                        n_genes++;
                }
                
                /* 2. Any gene positively regulated by a
                 *    currently active gene is turned on.
                 */
                for( unsigned int j = 0; j < on.size(); j++ ) {
                        on.at(j)->setConcentration( 1 );
                }
                
                /* 3. Any gene negatively regulated by a
                 *    currently active gene is turned off.
                 */
                for( unsigned int j = 0; j < off.size(); j++ ) {
                        off.at(j)->setConcentration( 0 );
                }       
        }

        data.close();
}

void executeNetworkT( AGraph* graph ) {
        vector<Specie*> on, off;
        vector<AVertex*> miRNA;
        AVertex *v, *vmiRNA;
        Specie* s1;
        Specie* s2;
        int n_species = 0, i;
        unsigned int j;
        string name, n;
        ofstream data1, data2, data3, data4, data5;
        bool stop, genestep = false;
        
        /*
        for( graph->it_first(); !graph->it_isDone(); graph->it_next() ) {
                n_species++;
        }
        */
        
        writeNetworkPlotSettings();
        
        ofstream file( (_plot_file+".input2").c_str() ); // ios::app 
        file << "set xrange[0:" << _num_steps+1 << "]\n";
        file << "set yrange[0:" << idgen.getCurrentID() << "]\n";
        file << "plot '" << _plot_file << "_gene.dat' with " << _plot_point_type << " lt -1";
        if( _plot_point_pt > 0 )
                file << " pt " << _plot_point_pt << "\n";
        file << "replot '" << _plot_file << "_protein.dat' with " << _plot_point_type << " lt 2";
        if( _plot_point_pt > 0 )
                file << " pt " << _plot_point_pt << "\n";
        file << "replot '" << _plot_file << "_mRNA.dat' with " << _plot_point_type << " lt 6";
        if( _plot_point_pt > 0 )
                file << " pt " << _plot_point_pt << "\n";
        file << "replot '" << _plot_file << "_miRNA.dat' with " << _plot_point_type << " lt 1";
        if( _plot_point_pt > 0 )
                file << " pt " << _plot_point_pt << "\n";
        file << "replot '" << _plot_file << "_ncRNA.dat' with " << _plot_point_type << " lt 3";
        if( _plot_point_pt > 0 )
                file << " pt " << _plot_point_pt << "\n";
        file.close();
        
        data1.open( (_plot_file+"_gene.dat").c_str() );
        data1 << "-1\t-1\n";
        data2.open( (_plot_file+"_protein.dat").c_str() );
        data2 << "-1\t-1\n";
        data3.open( (_plot_file+"_mRNA.dat").c_str() );
        data3 << "-1\t-1\n";
        data4.open( (_plot_file+"_miRNA.dat").c_str() );
        data4 << "-1\t-1\n";
        data5.open( (_plot_file+"_ncRNA.dat").c_str() );
        data5 << "-1\t-1\n";
        
        for( i = 0; i < _num_steps; i++ ) {
                n_species = 0;
                on.clear();
                off.clear();
                
                for( graph->it_first(); !graph->it_isDone(); graph->it_next() ) {
                        v = graph->it_current();
                        s1 = v->getData();
                        
                        //cout << s1->getName() << ":\n";
                        
                        if( s1->getConcentration() == 1 ) {
                                name = s1->getName();
                                //cout << name;
                                if( name.find("Gene") != string::npos ) {
                                        genestep = true;
                                        
                                        data1 << i << "\t" << s1->getID() << "\n";
                                        
                                        for( v->it_first(); !v->it_isDone(); v->it_next() ) {
                                                s2 = v->it_current();
                                                on.push_back( s2 );
                                                
                                                stop = false;
                                                for( j = 0; j < miRNA.size(); j++ ) {
                                                        vmiRNA = miRNA.at(j);
                                                        
                                                        for( vmiRNA->it_first(); !vmiRNA->it_isDone(); vmiRNA->it_next() ) {
                                                                if( vmiRNA->it_current() == s2 ) {      
                                                                        off.push_back( s2 );
                                                                        stop = true;
                                                                        break;
                                                                }
                                                        }
                                                        if( stop )
                                                                break;
                                                }
                                        }
                                } else if( name.find("Protein") != string::npos ) {
                                        data2 << i << "\t" << s1->getID() << "\n";
                                
                                        for( v->it_first(); !v->it_isDone(); v->it_next() ) {
                                                s2 = v->it_current();
                                
                                                if( v->currentWeight() > 0 ) {
                                                        on.push_back( s2 );
                                                        //cout << "\tactivates " << s2->getName() << endl;
                                                } else { 
                                                        off.push_back( s2 );
                                                        //cout << "\tstops " << s2->getName() << endl;
                                                }
                                        }
                                }
                                else if( name.find("mRNA") != string::npos ) {
                                        data3 << i << "\t" << s1->getID() << "\n";
                                        
                                        for( v->it_first(); !v->it_isDone(); v->it_next() ) {
                                                s2 = v->it_current();
                                                on.push_back( s2 );
                                        }
                                } else if( name.find("ncRNA") != string::npos ) {
                                        data5 << i << "\t" << s1->getID() << "\n";
                                        
                                        for( v->it_first(); !v->it_isDone(); v->it_next() ) {
                                                s2 = v->it_current();
                                                on.push_back( s2 );
                                        }
                                } else if( name.find("miRNA") != string::npos ) {
                                        data4 << i << "\t" << s1->getID() << "\n";
                                        miRNA.push_back( v );
                                
                                }
                                
                                /* 1. Any gene not regulated by a currently
                                 * active gene is turned off.
                                 */
                                s1->setConcentration( 0 );
                        }       
                        n_species++;
                }
                
                if( genestep ) {
                        miRNA.clear();
                        genestep = false;
                }
                
                /* 2. Any gene positively regulated by a
                 *    currently active gene is turned on.
                 */
                for( j = 0; j < on.size(); j++ ) {
                        on.at(j)->setConcentration( 1 );
                }
                
                /* 3. Any gene negatively regulated by a
                 *    currently active gene is turned off.
                 */
                for( j = 0; j < off.size(); j++ ) {
                        off.at(j)->setConcentration( 0 );
                }       
        }

        data1.close();
        data2.close();
        data3.close();
        data4.close();
        data5.close();
}

void initializeNetwork( const vector<Gene*>* genes, AGraph* graph ) {
        int r1 = -1, r3;
        float r2 = 0;
        unsigned int i;
        set<int> used_indexes;
        
        /* a particular gene is activated */
        if( _active_genes_type == 'g' ) {
                r1 = (int)_active_genes;
        } else {
                if( _active_genes >= 0 ) {
                        if( _active_genes > 1 ) {
                                cout << "Could not activate as many genes.\n";  
                                r2 = 1;
                        } else 
                                r2 = _active_genes;
                /* a gene with a high connectivity is activated (not necessarilly
                 * the one with the highest connectivity, just one with higher
                 * conectivity than the average) */
                } else if( _active_genes == -1 ) {
                        for( i = 0; i < genes->size(); i++ ) {
                                r2 = graph->find( (Specie *)genes->at(i) )->getConnectivity();
                                if( r2 >= graph->getAverageCon() ) {
                                        r1 = i;
                                        break;
                                }
                        }
                }
        }
        
        if( r1 >= 0 ) {
                genes->at(r1)->setConcentration( 1 );
        } else {
                /* random genes are activated */
                r3 = (int) (r2 * genes->size());
                for( int j = 0; j < r3; j++ ) {
                        /* Normal distribution */
                        if( _active_genes_normal_dist_mean >= 0 ) {
                                r1 = (int) mtrand.randNorm( _active_genes_normal_dist_mean * genes->size(), 
                                                                                        _active_genes_normal_dist_var * genes->size() );
                                if( r1 < 0 )
                                        r1 = 0;
                                else if( r1 >= genes->size() )
                                        r1 = genes->size() - 1;
                        /* Uniform distribution */
                        } else
                                r1 = mtrand.randInt(genes->size()-1); 
                        if( used_indexes.find( r1 ) == used_indexes.end() ) { 
                                genes->at(r1)->setConcentration( 1 );
                                used_indexes.insert( r1 );
                        } else
                                j--;
                }
        }
        
        used_indexes.clear();
}

AGraph* extractNetwork( AGraph* graph, const vector<Gene*>* genes, const vector<Protein*>* proteins, 
                                                const vector<Specie*>* mRNA, const vector<Specie*>* miRNA ) {
        unsigned int i, j, start_index, k, reg_size;
        int signal, positive_con = 0, negative_con = 0;
        int connectivity[proteins->size()];
        float bs = 0;
        string reg_region, is;
        string::size_type loc = string::npos;
        Protein* p;
        Gene* g;
        map<int,int> freq_input_aux, freq_output_aux, freq_input, freq_output;
        map<int,int> freq_mi_output_aux, freq_mi_output, freq_m_input_aux, freq_m_input;
        map<int,int>::iterator iter, iter_aux, end;
        
        writeStatisticsPlotSettingsToFile( "gene"+_input_con_file, "Number of inputs", "Number of genes" );
        writeStatisticsPlotSettingsToFile( "protein"+_output_con_file, "Number of outputs", "Number of proteins" );
        writeStatisticsPlotSettingsToFile( "miRNA"+_output_con_file, "Number of outputs", "Number of miRNAs" );
        writeStatisticsPlotSettingsToFile( "mRNA"+_input_con_file, "Number of inputs", "Number of mRNAs" );
        
        if( _algorithm != ARTGEN )      // TODO: sitio melhor
                initialize_binding_table();
         
        /* For each protein... */
        for( i = 0; i < proteins->size(); i++ ) {
                cout << "\rExtracting network... " << (i*100)/proteins->size() << "%" << flush;
                //cout << "." << flush;
                p = proteins->at(i);
                connectivity[i] = 0;
                freq_output_aux[i] = 0;
                
                if( p->getSequence().size() >= _binding_site_size ) {
                        /* ... search all genes for a match. */
                        for( j = 0; j < genes->size(); j++ ) {
                                g = genes->at(j);
                                reg_size = g->getRegRegionSize();
                                
                                if( reg_size < _binding_site_size )
                                        continue;
                                
                                iter = freq_input_aux.find(j);
                                if( iter == freq_input_aux.end() )
                                        freq_input_aux[j] = 0;
                                
                                if( _binding_site == 0 || reg_size < _binding_site_size )
                                        reg_region = g->getRegRegion();
                                else
                                        reg_region = g->getRegRegion().substr( reg_size - _binding_site_size, reg_size-1 );
                                start_index = 0;
                                
                                while( start_index <= p->getSequence().size() - _binding_site_size ) {  
                                        if( _algorithm == ARTGEN ) {
                                                loc = reg_region.find( p->getSequence().substr( start_index, _binding_site_size), 0 );
                                                //cout << "comparing " << p->getSequence().substr( start_index, _binding_site_size) <<
                                                //      " with reg region " << reg_region << endl;
                                        } else { // if( _algorithm == ONE_LEVEL )
                                                for( k = 0; k <= reg_size - _binding_site_size; k++ ) {
                                                        //cout << "k <= " << reg_size << " - " << _binding_site_size << " = " << reg_size - _binding_site_size << endl;
                                                        //cout << "comparing " << p->getSequence().substr( start_index, _binding_site_size) <<
                                                        //" with reg region " << reg_region.substr( k, _binding_site_size ) << endl;
                                                        //bs = dnaProteinBinding( p->getSequence().substr(start_index, _binding_site_size),
                                                        //                                          reg_region.substr( k, _binding_site_size ) );
                                                        bs = dnaProteinBinding( p->getSequence(), start_index, _binding_site_size, reg_region, k );
                                                        if( bs >= _binding_threshold )
                                                                break;
                                                        /*
                                                        else
                                                                cout << "Not found for " << p->getID() << " " 
                                                                << p->getSequence().substr(start_index, _binding_site_size) 
                                                                << " at " <<start_index
                                                                << " with size " << p->getSequence().size() 
                                                                << " on gene " << g->getID() << " " 
                                                                << reg_region.substr( k, _binding_site_size )
                                                                << " at " 
                                                                << k << " with strength "
                                                                << bs << endl;
                                                        */
                                                } 
                                        }
                                        
                                        if( loc != string::npos || bs >= _binding_threshold ) {
                                                freq_output_aux[i]++;
                                                freq_input_aux[j]++;
                                                /*
                                                cout << "Found match " 
                                                     << p->getSequence().substr(start_index, _binding_site_size) 
                                                     << " at position " << start_index << " in the protein " 
                                                     << p->getID() << " on gene " << g->getID() 
                                                     << " at position " << k << " with strength "
                                                     << bs << endl;
                                                */
                                                if( _inhibition_type == "random" )
                                                        signal = mtrand.rand() <= _inhibition_rate ? -1: 1;
                                                else
                                                        signal = p->getSignal( _inhibition );
                                                
                                                if( _algorithm == ARTGEN || _algorithm == ONE_LEVEL )
                                                        graph->AddEdge( p->getParent(), g, signal );
                                                else 
                                                        graph->AddEdge( p, g, signal );
                                                        
                                                if( signal > 0 )
                                                        positive_con++;
                                                else if( signal < 0 )
                                                        negative_con++;
                                                
                                                connectivity[i]++;      
                                                break;
                                        }
                                        start_index++; 
                                }
                        }
                }
        }
        
        /* Search for miRNA matches */
        if( _algorithm == TWO_LEVELS ) {
                for( i = 0; i < miRNA->size(); i++ ) {
                        freq_mi_output_aux[i] = 0;
                        
                        is = miRNA->at(i)->getSequence();
                        for( j = 0; j < mRNA->size(); j++ ) {
                                if( mRNA->at(j) != miRNA->at(i)->getParent() ) {
                                        loc = mRNA->at(j)->getSequence().find( is );
                                        if( loc != string::npos ) {
                                                iter = freq_m_input_aux.find(j);
                                                if( iter == freq_m_input_aux.end() )
                                                        freq_m_input_aux[j] = 0;
                                                freq_m_input_aux[j]++;
                                                
                                                graph->AddEdge( miRNA->at(i), mRNA->at(j), -1 );
                                                negative_con++;
                                                freq_mi_output_aux[i]++;
                                        }
                                }
                        }
                }
        }       
        
        for( iter_aux = freq_input_aux.begin(), end = freq_input_aux.end();  iter_aux != end;  iter_aux++ ) {
                //cout << "gene " << iter_aux->second
                iter = freq_input.find(iter_aux->second);
                if( iter == freq_input.end() )
                        freq_input[iter_aux->second] = 0;
                freq_input[iter_aux->second]++;
        }
        
        for( iter_aux = freq_output_aux.begin(), end = freq_output_aux.end();  iter_aux != end;  iter_aux++ ) {
                iter = freq_output.find(iter_aux->second);
                if( iter == freq_output.end() )
                        freq_output[iter_aux->second] = 0;
                freq_output[iter_aux->second]++;
        }
        
        for( iter_aux = freq_mi_output_aux.begin(), end = freq_mi_output_aux.end();  iter_aux != end;  iter_aux++ ) {
                iter = freq_mi_output.find(iter_aux->second);
                if( iter == freq_mi_output.end() )
                        freq_mi_output[iter_aux->second] = 0;
                freq_mi_output[iter_aux->second]++;
        }
        
        for( iter_aux = freq_m_input_aux.begin(), end = freq_m_input_aux.end();  iter_aux != end;  iter_aux++ ) {
                iter = freq_m_input.find(iter_aux->second);
                if( iter == freq_m_input.end() )
                        freq_m_input[iter_aux->second] = 0;
                freq_m_input[iter_aux->second]++;
        }
        
        ofstream data_input( ("gene"+_input_con_file+".dat").c_str() ); 
        for( iter= freq_input.begin(), end = freq_input.end();  iter != end;  iter++ ) 
        data_input << iter->first << "\t" << iter->second << endl;
        data_input.close();

        ofstream data_output( ("protein"+_output_con_file+".dat").c_str() ); 
        for( iter = freq_output.begin(), end = freq_output.end();  iter != end;  iter++ )
        data_output << iter->first << "\t" << iter->second << endl;
        data_output.close();
        
        data_output.open( ("miRNA"+_output_con_file+".dat").c_str() );
        for( iter = freq_mi_output.begin(), end = freq_mi_output.end();  iter != end;  iter++ )
        data_output << iter->first << "\t" << iter->second << endl;
        data_output.close();
        
        data_input.open( ("mRNA"+_input_con_file+".dat").c_str() ); 
        for( iter= freq_m_input.begin(), end = freq_m_input.end();  iter != end;  iter++ ) 
        data_input << iter->first << "\t" << iter->second << endl;
        data_input.close();
        
        graph->setPosCon( positive_con );
        //average(connectivity, proteins->size()) );
        graph->setNegCon( negative_con );
        return graph;  
}

AGraph* addVertexes( AGraph* graph, vector<Specie*>* species, string name ) {
        Specie* p;
        int count = 0;
        
        for( unsigned i = 0; i < species->size(); i++ ) {
                cout << "\rAdding the " << name << " to the graph... " << (i*100)/species->size() << "%" << flush;
                p = species->at(i);
                graph->AddVertex( p );
                if( p->getParent() != NULL ) {
                        //graph->AddVertex( p->getParent() );
                        graph->AddEdge( p->getParent(), p, 0 );
                }
                count++;
        }
        
        return graph;
}

vector<string::size_type> findStem( string& s ) {
        unsigned int i, j, k, l, count = 0;
        vector<string::size_type> pos;
        
        //cout << "Looking for stem in " << s << endl;
        for( i = _miRNA_binding_site_size - 1; i < s.size() - _miRNA_binding_site_size - 1; i++ ) {
                for( j = i; j < s.size(); j++ ) {
                        if( ( s[i] % 2 == 0 ? s[i]+1 : s[i]-1 ) == s[j] ) { 
                                k = i - 1;
                                l = j + 1;
                                count = 1;
                                while( k > 0 && ( s[k] % 2 == 0 ? s[k]+1 : s[k]-1 ) == s[l] ) {
                                        count++;
                                        k--;
                                        l++;
                                }
                                if( count >= _miRNA_binding_site_size ) {
                                        pos.push_back( k+1 );
                                        //cout << "found at " << k+1 << " - " << l-1 << endl;
                                        i = l-1;
                                        break;
                                }
                        }
                }
        }
        
        return pos;
}

vector<Specie*>* create_miRNA( vector<Specie*>* ncRNA ) {
        vector<Specie*>* miRNA = new vector<Specie*>;
        vector<string::size_type> pos;
        Specie* specie;
        string seq, mseq, id;
        unsigned int i, j;
        int cm = 0, cn = 0;
        map<int,int> freq_output_aux, freq_output;
        map<int,int>::iterator iter, iter_aux, end;
        
        writeStatisticsPlotSettingsToFile( "ncRNA"+_output_con_file, "Number of outputs", "Number of ncRNAs" );
        
        for( i = 0; i < ncRNA->size(); i++ ) {
                freq_output_aux[i] = 0;
                seq = ncRNA->at(i)->getSequence();
                
                if( seq.size() >= 2*_miRNA_binding_site_size ) {
                        pos = findStem( seq );
                        if( pos.size() > 0 ) {
                                for( j = 0; j < pos.size(); j++ ) {
                                        mseq = seq.substr( pos.at(j), _miRNA_binding_site_size );
                                        /*
                                        if( recycled ) {
                                                specie = new Specie();
                                                id = idgen.getCurrentID();
                                                specie->setID( id );
                                                //specie->setName( "miRNA" + currentID );
                                                specie->setName( "miRNA" + id ); 
                                                specie->setParent( ncRNA->at(i)->getParent() );
                                                specie->setSequence( mseq );
                                                miRNA->push_back( specie );
                                        } else {
                                                ncRNA->at(i)->setSequence( mseq );
                                                ncRNA->at(i)->setName( "miRNA" + ncRNA->at(i)->getID() ); 
                                        }
                                        */
                                        specie = new Specie();
                                        id = idgen.getCurrentID();
                                        specie->setID( id );
                                        specie->setName( "miRNA" + id );
                                        specie->setParent( ncRNA->at(i) );
                                        specie->setSequence( mseq );
                                        miRNA->push_back( specie );
                                        if( specie->getParent()->getName().find( "mRNA" ) != string::npos )
                                                cm++;
                                        else
                                                cn++;
                                        freq_output_aux[i]++;
                                }
                                //ncRNA->erase( ncRNA->begin() + i, ncRNA->begin() + i + 1 );
                                //i--;
                        }
                }
        }
        
        for( iter_aux = freq_output_aux.begin(), end = freq_output_aux.end();  iter_aux != end;  iter_aux++ ) {
                iter = freq_output.find(iter_aux->second);
                if( iter == freq_output.end() )
                        freq_output[iter_aux->second] = 0;
                freq_output[iter_aux->second]++;
        }
        
        ofstream data_output( ("ncRNA"+_output_con_file+".dat").c_str() ); 
        for( iter= freq_output.begin(), end = freq_output.end();  iter != end;  iter++ ) 
        data_output << iter->first << "\t" << iter->second << endl;
        data_output.close();
        
        /*
        for( i = 0; i < miRNA->size(); i++ )
                ncRNA->push_back( miRNA->at(i) );
        */
        cout << "Proportion of mRNA derived miRNAs: " << ((float)cm/(cm+cn))*100 << "\%" << endl;
        
        return miRNA;
}

string translate( string &seq ) {
        unsigned int i;
        int m, n, o;
        int c, add;
        vector<int>* v = new vector<int>;
        string p( "" );
        bool start_codon_found = false;
        
        if( seq.size() >= 3 ) {
                add = 1;

                for( i = 0; i < seq.size() - 2; i+=add ) {
                        //cout << "\tposition " << i << endl; 
                        m = seq[i] - 48;
                        n = seq[i+1] - 48;
                        o = seq[i+2] - 48;
                        //cout << "Codon " << m << n << o << endl;
                        c = translation_table[m][n][o];
                                
                        if( start_codon_found ) {
                                if( c == STOP )
                                        break;
                                v->push_back( c );
                        } else if( c == MET ) {
                                start_codon_found = true;
                                v->push_back( c );
                                add = 3;
                        }
                        //sprintf( buf, "%d", c );
                        //p += buf;
                }
                
                p = toProteinEncoding( v );
        }
        
        delete v;
        return p;
}

vector<Protein*>* translation( vector<Specie*>* genes, vector<Specie*>* ncRNA ) {
        vector<Protein*>* proteins = new vector<Protein*>();
        Protein* prot;
        Specie* g;
        string p, d;
        int msize, ncSize;
        map<int,int> freq, freq_nc, freq_m, freq_reg;
        map<int,int>::iterator iter, end;
        
        writeStatisticsPlotSettingsToFile( "mRNA"+_size_file, "mRNA size", "Number of mRNAs" );
        writeStatisticsPlotSettingsToFile( "protein"+_size_file, "Protein size", "Number of proteins" ) ;
        writeStatisticsPlotSettingsToFile( "ncRNA"+_size_file, "ncRNA size", "Number of ncRNAs" ) ;     
                
        /* ArtificialGenome translation:
         * p[i] = (g[i]+1) % 4
         */ 
        if( _algorithm == ARTGEN ) {
                char c[2];
                c[1] = '\n';
        
                for( unsigned int i = 0; i < genes->size(); i++ ) {
                        g = genes->at(i);
                        d = g->getSequence();
                        p = "";
                        for( unsigned int j = 0; j < d.size(); j++ ) {
                                c[0] = d[j]; 
                                c[0] = (( atoi(c) + 1 ) % _num_bases) + '0';
                                p += c[0];
                        }
                        prot = new Protein( idgen.getCurrentID(), g, p );
                        proteins->push_back( prot );
                }
        /* One additional level of regulation:
         * p[i] = genetic_code( g[i], g[i+1], g[i+2] )
         */
        } else if( _algorithm == ONE_LEVEL || _algorithm == TWO_LEVELS ) {
                initialize_translation_table();
                initialize_one_letter_protein_encoding_table();
                
                //cout << "\n\ttranslating " << genes->size() << " mRNAs\n";
                for( unsigned int i = 0; i < genes->size(); i++ ) {
                        msize = genes->at(i)->getSequence().size();
                        iter = freq_m.find(msize);
                        if( iter == freq_m.end() )
                                freq_m[msize] = 0;
                        freq_m[msize]++;
                        
                        g = genes->at(i);
                        d = g->getSequence();
                        p = translate(d);
                        
                        iter = freq.find(p.size());
                        if( iter == freq.end() )
                                freq[p.size()] = 0;
                        freq[p.size()]++;
                        
                        if( p.size() > 0 ) {
                                prot = new Protein( idgen.getCurrentID(), g, p );
                                proteins->push_back( prot );
                        } else { // it's a non-coding RNA
                                for( unsigned j = 0; j < ncRNA->size(); j++ ) {
                                        if( ncRNA->at(j)->getParent() == g ) {
                                                ncRNA->erase( ncRNA->begin() + j, ncRNA->begin() + j + 1 );
                                                j--;
                                        }
                                }
                                ncRNA->push_back( g );
                                genes->erase( genes->begin() + i, genes->begin() + i + 1 );
                                i--;
                        }
                }
        }
        /*
        cout << "\tyes = " << yes << " no = " << no << endl;
        cout << "\tgenerated " << proteins->size() << " proteins\n";
        cout << "\tmRNAs are now " << genes->size() << endl;
        cout << "\tncRNAs are now " << ncRNA->size() << endl;
        */
        
        for( unsigned j = 0; j < ncRNA->size(); j++ ) {
                ncSize = ncRNA->at(j)->getSequence().size();
                iter = freq_nc.find(ncSize);
                if( iter == freq_nc.end() )
                        freq_nc[ncSize] = 0;
                freq_nc[ncSize]++;
        }
        
        ofstream data( ("ncRNA"+_size_file+".dat").c_str() );
        for( iter = freq_nc.begin(), end = freq_nc.end();  iter != end;  iter++ )
        data << iter->first << "\t" << iter->second << endl;
        data.close();
        
        data.open( ("protein"+_size_file+".dat").c_str() ); 
        for( iter = freq.begin(), end = freq.end();  iter != end;  iter++ )
        data << iter->first << "\t" << iter->second << endl;
        data.close();
        
        data.open( ("mRNA"+_size_file+".dat").c_str() ); 
        for( iter = freq_m.begin(), end = freq_m.end();  iter != end;  iter++ )
        data << iter->first << "\t" << iter->second << endl;
        data.close();
        
        return proteins;
}

vector<Specie*>* splicing( const vector<Specie*>* mRNA ) {
        Search match;
        string::size_type iloc, floc;
        vector<Specie*>* ncRNA = new vector<Specie*>;
        Specie* s;
        string m, aux, temp, id;
        float f;
        map<int,int> freq2, freq_output_aux, freq_output;
        map<int,int>::iterator iter, iter_aux, end;
        
        writeStatisticsPlotSettingsToFile( "intron"+_size_file, "Introns size", "Number of introns" );
        writeStatisticsPlotSettingsToFile( "gene"+_output_con_file, "Number of outputs", "Number of genes" );
        
        f = (1-_U1_match) * (_left_U1.size() + _right_U1.size());
        match.setLimit( (int) (f+0.5) );

        for( unsigned int i = 0; i < mRNA->size(); i++ ) {
                freq_output_aux[i] = 0;
                m = mRNA->at(i)->getSequence();
                
                while(true) {
                        //cout << "Looking in mRNA " << i << ": " << m << endl;
                        iloc = match.find( m, _left_U1 );
                        if( iloc == string::npos )
                                break;
                        //cout << "\tfound left at " << iloc << endl;
                        aux = m.substr( iloc + _left_U1.size() );
                        floc = match.find( aux, _right_U1 );
                        if( floc == string::npos )
                                break;
                        floc += iloc + _left_U1.size();
                        //cout << "\tfound right at " << floc << endl;
                        
                        freq_output_aux[i]++;
                        s = new Specie();
                        id = idgen.getCurrentID();
                        s->setID( id );
                        s->setName( "ncRNA" + id );
                        s->setParent( mRNA->at(i)->getParent() );
                        s->setSequence( m.substr(iloc, floc-iloc) );
                        ncRNA->push_back( s );
                        aux = m.substr( 0, iloc );
                        //cout << "\taux " << aux << endl;
                        temp = m.substr( floc + _right_U1.size() );
                        //cout << "\ttemp " << temp << endl;
                        m.clear();
                        m = aux + temp;
                        
                        iter = freq2.find(s->getSequence().size());
                        if( iter == freq2.end() )
                                freq2[s->getSequence().size()] = 0;
                        freq2[s->getSequence().size()]++;
                }
                mRNA->at(i)->setSequence( m );
        }
        
        for( iter_aux = freq_output_aux.begin(), end = freq_output_aux.end();  iter_aux != end;  iter_aux++ ) {
                iter = freq_output.find(iter_aux->second);
                if( iter == freq_output.end() )
                        freq_output[iter_aux->second] = 0;
                freq_output[iter_aux->second]++;
        }
        
        ofstream data( ("intron"+_size_file+".dat").c_str() ); 
        for( iter = freq2.begin(), end = freq2.end();  iter != end;  iter++ )
        data << iter->first << "\t" << iter->second << endl;
        data.close();
        
        data.open( ("gene"+_output_con_file+".dat").c_str() ); 
        for( iter= freq_output.begin(), end = freq_output.end();  iter != end;  iter++ ) 
        data << iter->first << "\t" << iter->second << endl;
        data.close();
        
        return ncRNA;
}

vector<Specie*>* transcription( const vector<Gene*>* genes ) {
        vector<Specie*>* mRNA = new vector<Specie*>;
        Specie* sp;
        Gene* g;
        string t, seq, id;
        int b;
        unsigned int i, j;
        
        for( i = 0; i < genes->size(); i++ ) {
                t.clear();
                g = genes->at(i);
                sp = new Specie();
                id = idgen.getCurrentID();
                sp->setID( id );
                sp->setName( "mRNA" + id );
                sp->setParent( g );
                
                seq = g->getSequence();
                for( j = 0; j < seq.size(); j++ ) {
                        b = seq.at(j) - 48;
                        t += ( b % 2 == 0 ? b+1+48 : b-1+48 );
                }
                sp->setSequence( t );
                mRNA->push_back( sp );
        }
        
        return mRNA;
}

vector<Gene*>* createGenes( string &s ) {
        vector<Gene*>* genes = new vector<Gene*>();
        vector<string> tokens;
        vector<string::size_type> tok_loc;
        Gene* g;
        unsigned int i = 0, j, diff;
        unsigned int maxsize = 0;
        string t, test, sub;
        string::size_type prev_loc, loc, tloc, auxloc, add;
        char *aux, *token;
        char* delims = "][-*?!"; 
        char* v = new char[7];
        bool ok;
        Search match;
        float f;
        map<int,int> freq, freq_reg;
        map<int,int>::iterator iter, end;

        int count, count1;
        
        ofstream setfile( ("gene"+_size_file+".input").c_str() );
        setfile << "set title \"Gene size\"\n";
        setfile << "set xlabel \"Gene size\"\n";
        setfile << "set ylabel \"Number of genes\"\n";
        setfile << "set boxwidth 1\n";
        setfile << "set style fill solid border -1\n";
        setfile << "set terminal png\n";
        setfile << "set output '" << "gene"+_size_file+".png'\n";
        setfile << "load '" << "gene"+_size_file+".input2'\n";
        setfile.close();
        
        writeStatisticsPlotSettingsToFile( "reg"+_size_file, "Reg size", "Number of reg" );
                
        /* Search the genome for each of the defined promoters. */
        for( unsigned int k = 0; k < _promoter.size(); k++ ) {
                count1 = 0;
                aux = new char[_promoter.at(k).size()+1];
                
                /* Copy the promoter to an aux variable */
                strncpy( aux, _promoter.at(k).c_str(), _promoter.at(k).size()+1 );
                /* Search for short exact sequences in the promoter. */
                
                token = strtok( aux, delims );
                if( token == NULL ) 
                        cout << "createGenes: bad promoter\n";

                count = 0;
                while( token != NULL ) {
                        test = token;
                        tokens.push_back( test );
                        count += test.size();
                        token = strtok( NULL, delims );
                }
                free( aux );
                f = (1-_promoter_match) * count;
                //cout << "f is " << f << " and i is " << (int) (f+0.5) << endl;
                match.setLimit( (int) (f+0.5) );
                                        
                /* Find the positions on the promoter where those 
                 * sequences are. */
                for( i = 0; i < tokens.size(); i++ ) {
                        tok_loc.push_back( _promoter.at(k).find( tokens[i] ) );
                }
                //cout << "Token in position " << tok_loc << " in the promoter\n";
                
                /* Find the position on the genome where the first
                 * sequence is. */
                prev_loc = tok_loc[0];
                sub = s.substr(prev_loc);
                loc = match.find( sub, tokens[0] );
                add = prev_loc;
                
                i = 0;
                while( loc != string::npos ) {
                        loc += add;
                        //cout << "Token " << tokens[0] << " in position " << loc << " in the genome\n";
                        ok = true;
                        count1++;
                        
                        cout << "\rSearching the genome for genes... "
                             << (prev_loc*100)/_genome_size << "%" << flush;
                        
                        /* See if the spacing is ok */
                        for( j = 1; j < tokens.size(); j++ ) {
                                diff = tok_loc[j] - tok_loc[j-1];
                                //cout << "\tsearching for token " << tokens[j] << " between " 
                                //       << loc+diff << " and " << loc+diff+1 << endl;
                                if( loc+diff >= s.size() ) {
                                        ok = false;
                                        break;
                                }
                                sub = s.substr(loc+diff, tokens[j].size());
                                auxloc = match.find( sub, tokens[j] );

                                if( auxloc == string::npos ) {
                                        ok = false;
                                        break;
                                }
                                /*
                                } else {
                                        //cout << "\tother token found at " << auxloc << endl;  
                                }*/
                        }
                        
                        if( ok ) {
                                test = s.substr( loc - tok_loc[0], _promoter.at(k).size() );
                                //cout << "Substring of genome tested: " << test << " at " << loc - tok_loc[0] << endl;
                                
                                /* Fix sized genes */
                                if( _termination_seq == "none" ) {
                                        t = s.substr( (loc-tok_loc[0]) + _promoter.at(k).size(), _gene_size);
                                        //cout << "corresponding gene " << t << endl;
                                        
                                        /* Ignore last gene if it's size is different
                                         * from GENE_SIZE. */
                                        if( t.size() != _gene_size )
                                                break;
                                        
                                        g = new Gene( idgen.getCurrentID(), test, t, s.substr(prev_loc, (loc-tok_loc[0]) - prev_loc) );
                                        genes->push_back( g );
                                        
                                        if( t.size() > maxsize )
                                                maxsize = t.size();
                                        //data << i << "\t" << t.size() << endl;

                                        iter = freq.find(t.size());
                                        if( iter == freq.end() )
                                                freq[t.size()] = 0;
                                        freq[t.size()]++;
                                        
                                        iter = freq_reg.find(g->getRegRegionSize());
                                        if( iter == freq_reg.end() )
                                                freq_reg[g->getRegRegionSize()] = 0;
                                        freq_reg[g->getRegRegionSize()]++;
                                        
                                        prev_loc = (loc-tok_loc[0]) + _promoter.at(k).size() + _gene_size;
                                        //cout << "\tcontinuing search from " << prev_loc << " to " << s.size()-tokens[0].size()+1 << endl;
                                        sub = s.substr(prev_loc);
                                        loc = match.find( sub, tokens[0] );
                                        add = prev_loc;
                                        i++;
                                } else {
                                        /* Search for a termination sequence after the promoter. */
                                        tloc = s.find( _termination_seq, (loc-tok_loc[0]) + _promoter.at(k).size() );
                                        if( tloc != string::npos ) {
                                                t = s.substr( (loc-tok_loc[0]) + _promoter.at(k).size(), 
                                                                        tloc - ((loc-tok_loc[0]) + _promoter.at(k).size()) );
                                                
                                                g = new Gene( idgen.getCurrentID(), test, t, s.substr(prev_loc, (loc-tok_loc[0]) - prev_loc) );
                                                genes->push_back( g );
                                                
                                                if( t.size() > maxsize )
                                                        maxsize = t.size();
                                                
                                                //data << i << "\t" << t.size() << endl;
                                                
                                                iter = freq.find(t.size());
                                                if( iter == freq.end() )
                                                        freq[t.size()] = 0;
                                                freq[t.size()]++;
                                                
                                                iter = freq_reg.find(g->getRegRegionSize());
                                                if( iter == freq_reg.end() )
                                                        freq_reg[g->getRegRegionSize()] = 0;
                                                freq_reg[g->getRegRegionSize()]++;
                                                
                                                prev_loc = (loc-tok_loc[0]) + _promoter.at(k).size() +
                                                                   (tloc - ((loc-tok_loc[0]) + _promoter.at(k).size())) +
                                                                   _termination_seq.size();
                                                //cout << "\tcontinuing search from " << prev_loc << " to " << s.size()-tokens[0].size()+1 << endl;
                                                sub = s.substr(prev_loc);
                                                loc = match.find( sub, tokens[0] );
                                                add = prev_loc;
                                                i++;
                                        } else {
                                                sub = s.substr(loc+1);
                                                add = loc+1;
                                                loc = match.find( sub, tokens[0] );
                                        }
                                }
                        } else {
                                sub = s.substr(loc+1);
                                add = loc+1;
                                loc = match.find( sub, tokens[0] );
                        }
                        test.clear();
                        t.clear();
                }
                test.clear();
                tokens.clear();
                tok_loc.clear();
                //cout << "\n\tfound " << count1 << " candidates\n";
        }
        free( v );
        
        ofstream data( ("gene"+_size_file+".dat").c_str() );
        for( iter = freq.begin(), end = freq.end();  iter != end;  iter++ )
        data << iter->first << "\t" << iter->second << endl;
        data.close();
        
        ofstream file( ("gene"+_size_file+".input2").c_str() ); // ios::app 
        file << "set xrange[-1:" << maxsize << "]\n";
        file << "set yrange[0:]\n";
        file << "plot '" << "gene"+_size_file << ".dat' with boxes";
        file.close();

        data.open( ("reg"+_size_file+".dat").c_str() );
        for( iter = freq_reg.begin(), end = freq_reg.end();  iter != end;  iter++ )
        data << iter->first << "\t" << iter->second << endl;
        data.close();

        return genes;
}

AGraph* readNetwork( string& file ) {
        string line, word1, word2, word3, name;
    ifstream is;
    AGraph* graph = new AGraph;
    vector<Gene*>* genes = new vector<Gene*>;
    //vector<Specie*>*  species = new  vector<Specie*>;
    map<int,Specie*> species;
    
        int signal = 1;
        int max = 0;

    is.open(file.c_str());
    if( !is )
        cout << "Error reading network file." << endl;
        
        while( getline(is, line) ) {
        istringstream in(line);
        in >> word1;
        in >> word2;
        
        if(  word2.find( "label" ) != string::npos ) {
                if( atoi( word1.c_str() ) > max )
                        max = atoi( word1.c_str() );
        
                if( word2.find( "miRNA" ) != string::npos )
                        name = "miRNA"+word1;
                else if( word2.find( "mRNA" ) != string::npos )
                        name = "mRNA"+word1;
                else if( word2.find( "ncRNA" ) != string::npos )
                        name = "ncRNA"+word1;
                
                species[atoi( word1.c_str() )] = new Specie( word1, name );
                
                if( word2.find( "Protein" ) != string::npos )
                        species[atoi( word1.c_str() )] = new Protein( word1, NULL, "" );
                else if( word2.find( "Gene" ) != string::npos ) {
                        genes->push_back( new Gene( word1, "", "", "" ) );
                        species[atoi( word1.c_str() )] = genes->back();
                }
                graph->AddVertex( species[atoi( word1.c_str() )] );
        } else if( word2 == "[color=red]" ) {
                signal = -1;
        } else if( word2 == "[color=black]" ) {
                signal = 1;
        } else if( word2.find( "color" ) != string::npos ) {
                signal = 0;
        } else if( word2 == "->" ) {
                in >> word3;
                graph->AddEdge( species[atoi(word1.c_str())],
                                                species[atoi(word3.c_str())], signal );
        }
        }
        idgen.startAt( max+1 );
        initializeNetwork( genes, graph );
        
        is.close();     
        return graph;
}

/*
void printGeneVector( vector<Gene*>* g, ofstream &stat ) {
        stat << "** Printing gene vector *****************\n";
        for( unsigned int i = 0; i < g->size(); i++ ) 
                stat << *(g->at(i));
        stat << flush;
}

void printProteinVector( vector<Protein*>* g, ofstream &stat ) {
        stat << "** Printing protein vector *****************\n";
        for( unsigned int i = 0; i < g->size(); i++ )
                stat << *(g->at(i));
        stat << flush;
}
*/

void printSpecieVector( string n, vector<Specie*>* g, ofstream &stat ) {
        stat << "** Printing " << n << " vector *****************\n";
        for( unsigned int i = 0; i < g->size(); i++ )
                stat << *(g->at(i));
        stat << flush;
}

void printGenome( string& s, ofstream &stat ) {
        stat << "** Printing genome *****************\n";
        for( int i = 0; i < _genome_size; i++ ) {
                //cout << setw(3) << s[i];
                stat << s[i];
                if( i != 0 && i % 80 == 0 )
                        stat << endl;
        }
        stat << endl;
}

void printConfig( ostream& oss ) {
        oss << "Used configuration:" << endl;
        oss << "\tgenome size = " << _genome_size << "\talgorithm = " << _algorithm << endl;
        oss << "\tpromoter match = " << _promoter_match << "\tleft u1 = " << _left_U1 << endl;
        oss << "\tright u1 = " << _right_U1 << "\tu1 match = " << _U1_match << endl;
        oss << "\tgene size = " << _gene_size << "\tgenome_file = " << _genome_file << endl;
        oss << "\tnum bases = " << _num_bases << "\tbinding site size = " << _binding_site_size << endl;
        oss << "\tmiRNA bind site size = " << _miRNA_binding_site_size
                << "\tterm seq = " << _termination_seq << endl;
        oss << "\tnum steps = " << _num_steps << "\tactive genes = " << _active_genes_type 
                << " " << _active_genes << endl;
        oss << "\tnormal mean = " << _active_genes_normal_dist_mean << "\tnormal var = " 
                << _active_genes_normal_dist_var << endl;
        oss << "\tinhibition rate = " << _inhibition_rate << "\tbind threshold = " 
                << _binding_threshold << endl;
        oss << "\tbinding choice = " << _binding_choice << endl;
}

void printStatistics( ostream& oss, vector<Gene*>* genes, vector<Specie*>* mRNA, 
                                          vector<Specie*>* ncRNA, vector<Specie*>* miRNA, 
                                          vector<Protein*>* proteins, AGraph* graph ) {
        int poscon, negcon;
        
        oss << "Statistics: \n";
        oss << "\tfound " << genes->size() << " genes\n";
        if( mRNA != NULL )
                oss << "\tcreated " << mRNA->size() << " mRNAs\n";
        if( ncRNA != NULL )
                oss << "\tfound " << ncRNA->size() << " ncRNAs\n";
        if( miRNA != NULL )
                oss << "\tcreated " << miRNA->size() << " miRNAs\n";
        oss << "\tcreated " << proteins->size() << " proteins\n";
        poscon = graph->getPosCon();
        negcon = graph->getNegCon();
        oss.precision( 2 );
        oss << "\tpositive connections " << ((float)poscon/(poscon+negcon))*100 << "% (" << poscon << ")\n";
        oss << "\tnegative connections " << ((float)negcon/(poscon+negcon))*100 << "% (" << negcon << ")\n";
        oss << "\taverage connectivity " << graph->getAverageCon() << endl;
        oss << "\taverage gene open frame size " << endl;
        oss << "\taverage protein size " << endl;
}

void commandLineArgs( char* arg1, char* arg2 ) { 
        if( strcmp(arg1, "-geno_size") == 0 )
                _genome_size = atoi( arg2 );
        else if( strcmp(arg1, "-alg") == 0 ) {
                if( strcmp(arg2, "one_level") == 0 )
                        _algorithm = ONE_LEVEL;
                else if( strcmp(arg2, "two_levels") == 0) 
                        _algorithm = TWO_LEVELS;
        } else if( strcmp(arg1, "-prom") == 0 )
                _promoter.push_back( string( arg2 ) );
        else if( strcmp(arg1, "-prom_match") == 0 )
                _promoter_match = atof( arg2 );
        else if( strcmp(arg1, "-left_u1") == 0 )
                _left_U1 = string( arg2 );
        else if( strcmp(arg1, "-right_u1") == 0 )
                _right_U1 = string( arg2 );
        else if( strcmp(arg1, "-u1_match") == 0 )
                _U1_match = atof( arg2 );
        else if( strcmp(arg1, "-gene_size") == 0 )
                _gene_size = atoi( arg2 );
        else if( strcmp(arg1, "-geno_file") == 0 )
                _genome_file = string( arg2 );
        else if( strcmp(arg1, "-n_bases") == 0 )
                _num_bases = atoi( arg2 );
        else if( strcmp(arg1, "-bind_site_size") == 0 )
                _binding_site_size = atoi( arg2 );
        else if( strcmp(arg1, "-mirna_bind_site_size") == 0 )
                _miRNA_binding_site_size = atoi( arg2 );
        else if( strcmp(arg1, "-term_seq") == 0 )
                _termination_seq = string( arg2 );
        else if( strcmp(arg1, "-num_steps") == 0 )
                _num_steps = atoi( arg2 );
        else if( strcmp(arg1, "-act_genes") == 0 ) {
                _active_genes_type = 'r';
                _active_genes = atoi( arg2 );
        } else if( strcmp(arg1, "-norm_dist_mean") == 0 )
                _active_genes_normal_dist_mean = atof( arg2 );
        else if( strcmp(arg1, "-norm_dist_var") == 0 )
                _active_genes_normal_dist_var = atof( arg2 );
        else if( strcmp(arg1, "-inhib_rate") == 0 )
                _inhibition_rate = atof( arg2 );
        else if( strcmp(arg1, "-bind_thresh") == 0 )
                _binding_threshold = atof( arg2 );
        else if( strcmp(arg1, "-bind_choice") == 0 )
                _binding_choice = string( arg2 );
        else if( strcmp(arg1, "-verbose") == 0 ) {
                if( arg2[0] == 'n' )
                        cout.clear(ios::failbit);
        } else {
                cout << "Available options are:\n"
                         << "\t-geno_size -alg -prom -prom_match -left_u1 -right_u1\n"
                         << "\t-u1_match -gene_size -geno_file -n_bases -bind_site_size\n"
                         << "\t-mirna_bind_site_size -term_seq -num_steps -act_genes\n"
                         << "\t-norm_dist_mean -norm_dist_var -inhib_rate -bind_thresh\n"
                         << "\t-bind_choice -verbose" << endl;
                         exit(-1);
        }
}

void readConfig( char *f ) {
    string line, word;
    ifstream is;
        string prom;
        char verbose;

    is.open(f);
    if( !is )
        cout << "Error reading configuration file." << endl;
        
        _algorithm = ARTGEN;

    while( getline(is, line) ) {
        istringstream in(line);
        in >> word;
        if(word == "GENOME_SIZE") {
            in >> _genome_size;
        } else if(word == "ALGORITHM" ) {
                in >> prom;
                if( prom == "one_level" ) 
                        _algorithm = ONE_LEVEL;
                else if( prom == "two_levels" )
                        _algorithm = TWO_LEVELS;
                prom.clear();
        } else if(word == "PROMOTER") {
                in >> prom;
                _promoter.push_back( prom );
        } else if(word == "PROMOTER_MATCH") {
                in >> _promoter_match;
        } else if(word == "LEFT_U1") {
                in >> _left_U1;
        } else if(word == "RIGHT_U1") {
                in >> _right_U1;
        } else if(word == "U1_MATCH") {
                in >> _U1_match;
        } else if(word == "GENE_SIZE") {
                in >> _gene_size;
        } else if(word == "GENOME_FILE") {
                in >> _genome_file;
        } else if(word == "NUM_BASES") {
                in >> _num_bases;
        } else if(word == "BINDING_SITE_SIZE") {
                in >> _binding_site_size;
        } else if(word == "MIRNA_BINDING_SITE_SIZE" ) {
                in >> _miRNA_binding_site_size;
        } else if(word == "BINDING_SITE") {
                in >> _binding_site;
        } else if(word == "TERMINATION_SEQUENCE") {
                in >> _termination_seq;
        } else if(word == "NUM_STEPS") {
                in >> _num_steps;
        } else if(word == "GRAPH_FILE") {
                in >> _graph_file;
        } else if(word == "PLOT_FILE") {
                in >> _plot_file;
        } else if(word == "SIZE_FILE") {
                in >> _size_file;
        } else if(word == "INPUT_CON_FILE") {
                in >> _input_con_file;
        } else if(word == "OUTPUT_CON_FILE") {
                in >> _output_con_file;
        } else if(word == "ACTIVE_GENES") {
                in >> _active_genes_type;
                in >> _active_genes;
        } else if(word == "ACTIVE_GENES_NORMAL_DIST_MEAN") {
                in >> _active_genes_normal_dist_mean;
        } else if(word == "ACTIVE_GENES_NORMAL_DIST_VAR") {
                in >> _active_genes_normal_dist_var;
        } else if(word == "PLOT_POINT_TYPE") {
                in >> _plot_point_type;
        } else if(word == "PLOT_POINT_PT") {
                in >> _plot_point_pt;
        } else if(word == "INHIBITION_TYPE") {
                in >> _inhibition_type;
        } else if(word == "INHIBITION_RATE") {
                in >> _inhibition_rate;
        } else if(word == "INHIBITION") {
                in >> _inhibition;
        } else if(word == "BINDING_THRESHOLD") {
                in >> _binding_threshold;
        } else if(word == "BINDING_CHOICE") {
                in >> _binding_choice;
        } else if(word == "READ_NETWORK_FILE") {
                in >> _read_network_file;
        } else if(word == "SAVE_NETWORK_FILE") {
                in >> _save_network_file;
        } else if(word == "STAT_FILE" ) {
                in >> _stat_file;
        } else if(word == "VERBOSE" ) {
                in >> verbose;
                if( verbose != 'y' )
                        cout.clear(ios::failbit);
        } else if(word.at(0) != '#') {
                cout << "Unrecognized configuration line: " << word << endl;
        }
    }

    is.close();
}

int main(int argc, char *argv[]) {
        string gene_seq( "" );
        char c;
        vector<Gene*>* genes;
        vector<Protein*>* proteins = NULL;
        vector<Specie*>* mRNA = NULL;
        vector<Specie*>* ncRNA = NULL;
        vector<Specie*>* miRNA = NULL;
        AGraph* graph = new AGraph();
        ofstream stat;
        
        initialize_translation_table();
        
        /* Read the configuration file */
        if( argc < 2 ) {
                cout << "Missing configuration file." << endl;
                return 1;
        }
        readConfig( argv[1] );
        
        for( int i = 2; i < argc; i = i + 2 )
                commandLineArgs( argv[i], argv[i+1] );
        cout << "Reading configuration file... done\n";
        
        stat.open( _stat_file.c_str() );
        if( !stat ) 
                cout << "Could not open " << _stat_file << " for writing\n";
        printConfig( cout );
        printConfig( stat );
        
        if( _read_network_file != "none" ) {
                cout << "Reading network from file... " << flush;
        AGraph* graph = readNetwork( _read_network_file );
        cout << "done\n";
                cout << "Execute network... " << flush;
                executeNetworkT( graph );
                cout << "done\n";
                delete graph;
        return 0;
        }
        
        /* Generate the genome */
        if( _genome_file == "none" ) {
                cout << "Generating the genome... " << flush;
                for( int i = 0; i < _genome_size; i++ ) {
                        c = mtrand.randInt(_num_bases-1) + '0';
                        gene_seq += c;
                }
                cout << "done\n";
                printGenome( gene_seq, stat );  // print the genome only if it's random
        } else {
                cout << "Reading the genome from " << _genome_file << "... ";
                ifstream is;
                string line;
                is.open( _genome_file.c_str() );
                while( getline( is, line ) ) {
                        gene_seq += line;
                }
                _genome_size = gene_seq.size();
                is.close();
                cout << "done\n";
        }
        cout << "size: " << _genome_size << endl;
        
        cout << "Searching the genome for genes... " << flush;
        genes = createGenes( gene_seq );
        cout << "\rSearching the genome for genes... done\n";
        gene_seq.clear();
        
        if( genes->size() <= 0 ) {
                cout << "No genes found... quitting.\n";
                return 0;
        }
        //printGeneVector( genes, stat );
        printSpecieVector( "gene", (vector<Specie*>*)genes, stat );
        cout << "\tfound " << genes->size() << " genes\n";
        
        if( _algorithm == ARTGEN || _algorithm == ONE_LEVEL ) {
                cout << "Generating proteins... " << flush;
                proteins = translation( (vector<Specie*>*)genes, NULL );
                cout << "done\n";
                //printProteinVector( proteins, stat );
                printSpecieVector( "protein", (vector<Specie*>*)proteins, stat );
                
                cout << "Adding the genes to the graph... " << flush;
                graph = addVertexes( graph, (vector<Specie*>*)genes, "genes" );
                cout << "\rAdding the genes to the graph... done\n"; 
                cout << "Extracting network... " << flush;
                graph = extractNetwork( graph, genes, proteins, NULL, NULL );
                cout << "\rExtracting network... done\n";
                cout << "Printing graph to " << _graph_file << ".dot" << "... " << flush;
                graph->print( _graph_file, _algorithm );
                cout << "done\n";
                        
                cout << "Initialize network... " << flush;
                initializeNetwork( genes, graph );
                cout << "done\n";
                
                cout << "Execute network... " << flush;
                executeNetwork( graph );
                cout << "done\n";
        } else if( _algorithm == TWO_LEVELS ) {
                cout << "Generating primary transcript... " << flush;
                mRNA = transcription( genes );
                cout << "done\n";
                printSpecieVector( "pre-mRNA", mRNA, stat );
                
                cout << "Splicing primary transcript... " << flush;
                /* mRNA is modified in this function */
                ncRNA = splicing( mRNA );
                cout << "done\n";
                printSpecieVector( "mRNA", mRNA, stat );
                printSpecieVector( "ncRNA", ncRNA, stat );
                
                cout << "Generating proteins... " << flush;
                proteins = translation( mRNA, ncRNA );
                cout << "done\n";
                printSpecieVector( "updated mRNA", mRNA, stat );
                printSpecieVector( "updated ncRNA", ncRNA, stat );
                //printProteinVector( proteins, stat );
                printSpecieVector( "protein", (vector<Specie*>*)proteins, stat );       
                cout << "\tcreated " << proteins->size() << " proteins\n";
                
                cout << "Generating miRNAs... " << flush;
                miRNA = create_miRNA( ncRNA );
                cout << "done\n";
                cout << "\tcreated " << miRNA->size() << " miRNAs\n";
                printSpecieVector( "miRNAs", miRNA, stat );
                
                cout << "Adding the genes to the graph... " << flush;
                graph = addVertexes( graph, (vector<Specie*>*)genes, "genes" );
                cout << "\rAdding the genes to the graph... done\n";
                cout << "Adding the mRNAs to the graph... " << flush;
                graph = addVertexes( graph, mRNA, "mRNAs" );
                cout << "\rAdding the mRNAs to the graph... done\n";
                cout << "Adding the proteins to the graph... " << flush;
                graph = addVertexes( graph, (vector<Specie*>*)proteins, "proteins" );
                cout << "\rAdding the proteins to the graph... done\n";
                cout << "Adding the ncRNAs to the graph... " << flush;
                graph = addVertexes( graph, ncRNA, "ncRNAs" );
                cout << "\rAdding the ncRNAs to the graph... done\n";
                cout << "Adding the miRNAs to the graph... " << flush;
                graph = addVertexes( graph, miRNA, "miRNA" );
                cout << "\rAdding the miRNAs to the graph... done\n";
                
                cout << "Extracting network... " << flush;
                graph = extractNetwork( graph, genes, proteins, mRNA, miRNA );
                cout << "\rExtracting network... done\n";
                
                cout << "Printing graph to " << _graph_file << ".dot" << "... " << flush;
                graph->print( _graph_file, _algorithm );
                cout << "done\n";
                        
                cout << "Initialize network... " << flush;
                initializeNetwork( genes, graph );
                cout << "done\n";
                
                cout << "Execute network... " << flush;
                executeNetworkT( graph );
                cout << "done\n";
        }
        
        printStatistics( cout, genes, mRNA, ncRNA, miRNA, proteins, graph );
        printStatistics( stat, genes, mRNA, ncRNA, miRNA, proteins, graph );
        
        //delete graph;
        delete genes;
        delete proteins;
        if( mRNA != NULL )
                delete mRNA;
        if( ncRNA != NULL )
                delete ncRNA;
        
        stat.close();
        
        return 0;
}

