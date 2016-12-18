#include "MEMsReader.hpp"
#include "MEMsList.hpp"

MemsReader::MemsReader(const std::string& fpath) {
    this->fpath = fpath;
}

void MemsReader::addPattern(const std::string& pattern_id, const int& pattern_length, std::forward_list<std::vector<int> > MEMs) {
    unsigned int UI_MemoryUsed = GetMemoryTotalPhys();
    std::cout << UI_MemoryUsed << std::endl;
    MemsList ml (pattern_length);
    std::cout << "Adding pattern..." << std::endl;
    while(!MEMs.empty()) {
	std::vector<int> m = MEMs.front();	
	MEMs.pop_front();
	ml.addMem(m[0], m[1], m[2]);
    }
    /**
    for(std::forward_list<std::vector<int> >::iterator it=MEMs.begin(); it != MEMs.end(); ++it) {
	std::cout << "0" << std::endl;
	std::vector<int> m = (*it);
	ml.addMem(m[0], m[1], m[2]);
    }
    **/
    std::cout << "Making pair..." << std::endl;
    std::pair<std::string, MemsList&> p (pattern_id, ml);
    std::cout << "Pushing..." << std::endl;
    patterns.push_front(p);
}

std::vector<int> MemsReader::extractMEM(std::string line) {
    std::vector<int> m (3);
    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;
    int i = 0;
    while ((pos = line.find(delimiter)) != std::string::npos) {
	token = line.substr(0, pos);
	m[i] = stoi(token);
	i++;
	line.erase(0, pos + delimiter.length());
    }
    m[i] = stoi(line);
    return m;
}

void MemsReader::readMEMsFile() {
    patterns.clear();
    bool flag = false;
    std::string line;
    std::string pattern_id;
    int pattern_length;
    std::forward_list<std::vector<int> > MEMs;
    std::ifstream memsFile(fpath);
    int i = 1;
    if (memsFile.is_open()) {
	while(getline(memsFile,line)) {
	    if(line.length() != 0) {
		if(line.at(0) == '#') {
		    if(flag) {
			std::cout << ++i << std::endl;
			addPattern(pattern_id, pattern_length, MEMs);
			std::cout << "-------------------" << std::endl;
			MEMs.clear();
			flag = false;
		    }
		    else {
			if(line.substr(0,7).compare("# P.len") == 0) {
			    pattern_length = stoi(line.substr(line.length()-2,line.length()));
			}
		    }
		}
		if(flag) {
		    std::cout << line << std::endl;
		    std::vector<int> mem = extractMEM(line);
		    MEMs.push_front(mem);
		}
		if(line.at(0) == '>') {
		    pattern_id = line.substr(2, line.length());
		    std::cout << pattern_id << std::endl;
		    flag = true;
		}
	    }
	}
	memsFile.close();
    }
    else {
	std::cout << "Unable to open MEMs file" << std::endl;
    }
}

void MemsReader::print() {
    for(std::forward_list<std::pair<std::string, MemsList> >::iterator it=patterns.begin(); it != patterns.end(); ++it) {
	std::pair<std::string, MemsList> pattern = (*it);
	std::cout << pattern.first << std::endl;
	MemsList MEMs = pattern.second;
	int i = 0;
	while(i<MEMs.getLength()) {
	    std::forward_list<Mem> l = MEMs.getMems(i);
	    std::cout << i << std::endl;
	    for(std::forward_list<Mem>::iterator it1=l.begin(); it1 != l.end(); ++it1) {
		Mem m = (*it1);
		std::cout << m.toStr() << std::endl;
	    }
	    i++;
	}
    }
}
