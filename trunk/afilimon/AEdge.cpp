#include "AEdge.h"
#include "AVertex.h"
/* Class AVertex is just used as a pointer reference in
 * AEdge.h. Thus a forward declaration is sufficient. 
 * But  AEdge.cpp uses class AVertex in substance so it
 * explicitly includes AVertex.h.
 */

AEdge::AEdge(AVertex *Vert, AEdge *nextEdge, float w) {
        end = Vert;
        next = nextEdge;
        weight = w;
}

AEdge::~AEdge() {
        delete next;
}

AVertex *AEdge::getEnd() {
        return end;
}

void AEdge::print( AVertex* v, ofstream &myfile ) {     
        //cout << weight << end->getData()->getName() << " ";
  
        if( weight == 0 ) 
                myfile << "edge [color=gray]\n";
        else if( weight < 0 )
                myfile << "edge [color=red]\n";
        else
                myfile << "edge [color=black]\n";
        myfile << v->getData()->getID() << " -> ";
        myfile << end->getData()->getID() << "\n";
        // descomentar para imprimir o peso
        // << "[label=\"" << weight << "\"]\n";
  
        if( next != NULL )
        next->print( v, myfile );
}

AEdge *AEdge::getNext() {
        return next;
}

float AEdge::getWeight() {
        return weight;
}

