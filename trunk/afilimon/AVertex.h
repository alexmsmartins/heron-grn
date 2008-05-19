#ifndef AVERTEX_H_
#define AVERTEX_H_

#ifndef AEDGE_H
#include "AEdge.h"
#endif
#include "Specie.h"
#include <iostream>
#include <fstream>

class AVertex {
        public:
        AVertex( Specie*, AVertex* );
        ~AVertex();
        //string getData();
        Specie* getData();
        void printEdges(ofstream&, ofstream&, int);
        void printGraph(ofstream&, ofstream&, int);
        void connectTo(AVertex*, int);
        AVertex* getNext();
        AEdge* getFirstEdge();
        int getConnectivity();
        
        Specie* it_current();
        bool it_isDone();
        void it_first();
        void it_next();
        float currentWeight();
 private:
  Specie* data;
  //string data;
  AEdge* edges;
  AVertex* next;
  
  AEdge* current;
  int connectivity;
};

#endif /*AVERTEX_H_*/


