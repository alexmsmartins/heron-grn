#include "AVertex.h"

AVertex::AVertex(Specie* theData, AVertex* nextVert) {
  data = theData;
  next = nextVert;
  edges = NULL;     // no edges 
  connectivity = 0;
}

AVertex:: ~AVertex() {
        delete data;
        delete next;
        delete edges;
}

Specie* AVertex::getData() {
  return data;
}
  
void AVertex::printEdges( ofstream& myfile, ofstream& aux, int algorithm ) {
        string color;
        string name;
        
        name = getData()->getName();
        if( name.find("Protein") != string::npos )
                color = "green";
        else if( name.find("Gene") != string::npos )
                color = "black";
        else if( name.find("mRNA") != string::npos )
                color = "gray";
        else if( name.find("ncRNA") != string::npos )
                color = "blue";
        else if( name.find("miRNA") != string::npos )
                color = "yellow";
                        
        myfile << getData()->getID() << " [label=\"" << getData()->getName() <<"\",color=" 
                   << color << "]\n";
        
        if( edges != NULL ) 
                edges->print( this, aux );
}  

void AVertex::printGraph( ofstream& file, ofstream& aux, int algorithm ) {
  printEdges( file, aux, algorithm );
  if( next != NULL )
    next->printGraph( file, aux, algorithm );
}
  
void AVertex::connectTo(AVertex *Vert, int weight ) {
  // allocate memory for a new Edge, set its Vertex pointer to point 
  // to Vert, and its Edge pointer to point to the rest of edges 
  AEdge *newEdge = new AEdge(Vert, edges, weight);

  // make the new edge the first edge of the vertex
  edges = newEdge; 
  connectivity++;
}

AVertex *AVertex::getNext() {
  return next;
}

AEdge *AVertex::getFirstEdge() {
  return edges;
}

int AVertex::getConnectivity() {
        return connectivity;
}

void AVertex::it_first() {
        current = edges;
}

bool AVertex::it_isDone() {
        if( current == NULL )
                return true;
        return false;
}

void AVertex::it_next() {
        current = current->getNext();
}

Specie* AVertex::it_current() {
        return current->getEnd()->getData();
}

float AVertex::currentWeight() {
        return current->getWeight();
}

