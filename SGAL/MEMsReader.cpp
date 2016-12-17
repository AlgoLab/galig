#include "MEMsReader.hpp"
//#include "MEMsList.hpp"

MemsReader::MemsReader() { }

void MemsReader::addPattern(const std::string& pattern_id, const int& pattern_length, std::forward_list<std::vector<int> > MEMs) {
    MemsList ml (pattern_length);

    for(std::forward_list<std::vector<int> >::iterator it=MEMs.begin(); it != MEMs.end(); ++it) {
	std::vector<int> m = (*it);
	ml.addMem(m[0], m[1], m[2]);
    }
    std::pair<std::string, MemsList> p (pattern_id, ml);
    patterns.push_front(p);
}

void MemsReader::getPatterns() {
    std::pair<std::string, MemsList> p = patterns.front();
    MemsList ml = p.second;
    int i = 0;
    while(i<ml.getLength()) {
	std::forward_list<Mem> l = ml.getMems(i);
	std::cout << i << std::endl;
	for(std::forward_list<Mem>::iterator it=l.begin(); it != l.end(); ++it) {
	    Mem m = (*it);
	    std::cout << m.toStr() << std::endl;
	}
	i++;
    }
}
