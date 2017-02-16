#include "MEMsGraph.hpp"

MemsGraph::MemsGraph(SplicingGraph& graph, std::list<Mem>& MEMs, const int& L) {
}
/**  
  Graph = TNGraph::New();
  addNode(Mem(0,0,0));
  int curr_p = 1;
  plen = ml.getLength();
  this->perc = perc;
  this->K = K;
  while(curr_p < plen) {
	std::forward_list<Mem> mems1 = ml.getMems(curr_p);
	for(auto it1=mems1.begin(); it1!=mems1.end(); ++it1) {
      Mem m1 = (*it1);
      if(!isNode(m1)) {
		double u = (((100 - this->perc) * plen) / 100) + 1;
		if(m1.p <= u + K) {
          addNode(m1);
          addEdge(Mem(0,0,0), m1, 0);
		}
		else {
          continue;
		}
      }
      //int i = m1.p + 1;
      int i = m1.p + m1.l - K;
      if(i<=m1.p) {
		i = m1.p + 1;//m1.p + m1.l;
      }
      while(i < plen && i <= m1.p + m1.l + K) {
		std::forward_list<Mem> mems2 = ml.getMems(i);
		for(auto it2=mems2.begin(); it2!=mems2.end(); ++it2) {
          Mem m2 = (*it2);
          //std::cout << "Checking " << m1.toStr() << " -> "<< m2.toStr() << std::endl;
          if(m1.p + m1.l < m2.p + m2.l) {
			if(g.rank(m1.t - 1) == g.rank(m2.t - 1)) {
              //Stesso esone
              if(m2.t > m1.t && m2.t <= m1.t + m1.l + K && m1.t + m1.l < m2.t + m2.l) {
				if(!isNode(m2)) {
                  //std::cout << "1 Adding " << m2.toStr() << std::endl;
                  addNode(m2);
				}
				int wt = m2.t - m1.t - m1.l;
				int wp = m2.p - m1.p - m1.l;
				int w;
				if(wt<0 || wp<0) {
				  w = abs(wt - wp);
				}
				else {
				  w = max(wt, wp);
				}
				//std::cout << "1 Adding " << m1.toStr() << " -> " << m2.toStr() << std::endl;
				addEdge(m1, m2, w);
              }
			} else {
              //Esoni diversi
              //std::cout << g.rank(m1.t-1) << ", " << g.rank(m2.t-1) << std::endl;
              std::vector<int> curr_edge { g.rank(m1.t-1), g.rank(m2.t-1) };
              if(g.contain(curr_edge)) {
				//std::cout << m1.toStr() << " " << m2.toStr() << ", " << g.select(g.rank(m1.t-1) + 1) << " " << g.select(g.rank(m2.t-1)) << std::endl;
				if(m1.t + m1.l >= g.select(g.rank(m1.t-1) + 1) + 1 - K && m2.t <= g.select(g.rank(m2.t-1)) + 1 + 1 + K) {
                  if(!isNode(m2)) {
					//std::cout << "2 Adding " << m2.toStr() << std::endl;
					addNode(m2);
                  }
                  int wt = (g.select(g.rank(m1.t-1) + 1) - m1.t - m1.l) + (m2.t - g.select(g.rank(m2.t-1)) - 1);
                  int wp = abs(m2.p - m1.p - m1.l);
                  int w = max(wt, wp);
                  addEdge(m1, m2, w);
                  //std::cout << "2 Adding " << m1.toStr() << " -> " << m2.toStr() << std::endl;
				}
              }
			}
          }
		}
		i++;
      }
	}
	curr_p++;
  }
  subpaths = std::vector<std::vector<std::vector<int> > >(Graph->GetNodes(), { std::vector<std::vector<int> > { std::vector<int> { } } });
}
**/
