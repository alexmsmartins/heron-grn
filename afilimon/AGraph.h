#ifndef AGRAPH_H_
#define AGRAPH_H_

#include "AVertex.h"
#include "AEdge.h"

class AGraph {
public:
        AGraph();
        ~AGraph();
        AVertex* find( Specie* );
        bool AddVertex( Specie* );
        bool AddEdge( Specie*, Specie*, int );
        void print( string, int );
        
        void setPosCon( int );
        void setNegCon( int );
        int getPosCon();
        int getNegCon();
        float getAverageCon();
        int getNumVert();
        
        AVertex* it_current();
        void it_next();
        void it_first();
        bool it_isDone();
private: 
        AVertex *first;  // a pointer to the first vertex of the graph
                     // NULL if the graph does not have vertexes
    AVertex* current;
    int pos_connections;        // activation
    int neg_connections;        // inhibition
    int nvert;          // number of vertexes in the graph
};
#endif /*AGRAPH_H_*/

AGraph::AGraph() {
        first = NULL;
        pos_connections = 0;
        neg_connections = 0;
        nvert = 0;
}
 
AGraph::~AGraph() {
        delete first;
}

AVertex *AGraph::find(Specie* theData) {
        AVertex *vPtr = first;
        
        while (vPtr != NULL) {
                if ( vPtr->getData() == theData )
                  return vPtr;
                vPtr = vPtr->getNext();
        }
        
        return NULL;
}

bool AGraph::AddVertex(Specie* theData) { 
        if ( find(theData) != NULL )
        return false;

        // allocate memory for new vertex with data theData,
        // make it point to the previous first vertex
        AVertex *newVertex = new AVertex(theData, first);
        // make the new vertex the first one in the list of vertexes
        first = newVertex;
        nvert++;
        return true;
}

bool AGraph::AddEdge(Specie* Begin, Specie* End, int weight) {
        // find the pointer to the end vertex
        AVertex *vEnd = find(End);
        
        // if the vertex is not in the graph, cannot add the edge
        if (vEnd == NULL) return false;
        
        // find the pointer to the start vertex
        AVertex *vBegin = find(Begin);
        
        // if the vertex is not in the graph, cannot add the edge
        if (vBegin == NULL) return false;
        
        // connect the start vertex to the end one
        vBegin->connectTo(vEnd, weight);
        return true;
}
  
void AGraph::print( string file, int algorithm ) {
        ofstream myfile, auxfile;
        file += ".dot";
  
        if (first == NULL) 
        cout << "Graph has no vertexes " << endl;
        else {
                myfile.open ( file.c_str() );
                myfile << "digraph simple {\n\n";
        
                // clear the file
                auxfile.open ("aux.dot");
        
                first->printGraph( myfile, auxfile, algorithm );
                auxfile.close();
                ifstream myfile2( "aux.dot" );
                char ch;
                while( myfile2.get(ch) )
                        myfile.put(ch); 
                myfile << "\n}";
                myfile2.close();
                myfile.close();
                if( remove( "aux.dot" ) != 0 )
                cout << "Error deleting aux.dot\n";
        }
}

void AGraph::setPosCon( int a ) {
        pos_connections = a;
}

void AGraph::setNegCon( int a ) {
        neg_connections = a;
}

int AGraph::getPosCon() {
        return pos_connections;
}

int AGraph::getNegCon() {
        return neg_connections;
}

float AGraph::getAverageCon() {
        return (pos_connections + neg_connections) / nvert;
}

int AGraph::getNumVert() {
        return nvert;
}

bool AGraph::it_isDone() {
        if( current == NULL )
                return true;
        return false;
}

void AGraph::it_next() {
        current = current->getNext();
}

AVertex* AGraph::it_current() {
        return current;
}

void AGraph::it_first() {
        current = first;
}

