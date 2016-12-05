/////////////////////////////////////////////////
// SNAP library
#include "stdafx.h"
#include "Snap.h"

#include "base.cpp"
#include "gnuplot.cpp"
#include "linalg.cpp"

#include "gbase.cpp"
#include "util.cpp"
#include "attr.cpp"     		 // sparse attributes

// graph data structures
#include "graph.cpp"         // graphs
#include "graphmp.cpp"       // graphs
//#include "mmgraph.cpp"       // multimodal graphs
#include "network.cpp"       // networks
#include "networkmp.cpp"     // networks OMP
#include "timenet.cpp"       // time evolving networks
#include "mmnet.cpp"         // multimodal networks

// algorithms
#include "subgraph.cpp"      // subgraph manipulations
#include "anf.cpp"           // approximate diameter calculation
#include "cncom.cpp"         // connected components
#include "alg.cpp"           // misc graph algorithms
#include "gsvd.cpp"          // SVD and eigenvector computations
#include "gstat.cpp"         // graph statistics
#include "centr.cpp"         // centrality measures
#include "cmty.cpp"          // community detection algorithms
#include "flow.cpp"          // network flow algorithms
#include "coreper.cpp"       // core-periphery algorithms
#include "triad.cpp"         // clustering coefficient and triads
#include "casc.cpp"          // cascade detection^M

// graph generators
#include "ggen.cpp"          // graph generators (preferential attachment, small-world, ...)
#include "ff.cpp"            // forest fire graph generator

#include "gviz.cpp"
#include "ghash.cpp"
#include "statplot.cpp"
#include "gio.cpp"

// table data structures and algorithms
#include "table.cpp"         // table
#include "conv.cpp"
#include "numpy.cpp"         // numpy conversion

