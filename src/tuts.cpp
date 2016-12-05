#include <iostream>

#include "Snap.h"
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

int main() {
  cout << "Tuts" << endl;
  
  TPt<TNodeEDatNet<TInt, TInt> >  Net = TNodeEDatNet<TInt, TInt>::New();
  Net->AddNode(1,3);
  Net->AddNode(2,2);
  Net->AddNode(3,1);
  
  Net->AddEdge(1, 2, 5);
  Net->AddEdge(2, 3, 2);
  
  for (TNodeEDatNet<TInt, TInt>::TNodeI NI = Net->BegNI(); NI < Net->EndNI(); NI++) { 
    printf("node id %d with %d\n", NI.GetId(), NI.GetDat());
  }
  
  for (TNodeEDatNet<TInt, TInt>::TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++) { 
    printf("edge (%d, %d, %d)\n", EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat()); 
  }
  
  TSnap::DrawGViz(Net, gvlDot, "graph.png", "");

}
