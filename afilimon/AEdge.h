#ifndef AEDGE_H_
#define AEDGE_H_

#include <iostream>
#include <fstream>

using namespace std;

/* Class AVertex is included only as a pointer reference. 
 * The size or actual content of AVertex are not important to
 * AEdge.h or  AEdge.cpp. Thus only a forward declaration has 
 * been included in AEdge.h.
 */
class AVertex;

class AEdge{
 public:
  AEdge(AVertex*, AEdge*, float);
  ~AEdge();
  AVertex *getEnd();
  void print(AVertex*, ofstream&);
  AEdge *getNext();
  float getWeight();

 private:
  AVertex *end;
  AEdge *next;
  float weight;
};

#endif /*AEDGE_H_*/

