#include "MEMsList.hpp"

MemsList::MemsList(const int& l) {
    length = l+1;
    std::vector<std::forward_list<Mem> > mems (length, std::forward_list<Mem> ());
    this->mems = mems;
}

void MemsList::addMem(const Mem& mem) {
    mems[mem.p].push_front(mem);
}

std::forward_list<Mem> MemsList::getMems(const int& i) {
    return mems[i];
}

int MemsList::getLength() {
    return length - 1;
}

void MemsList::print() {
    int i = 0;
    while(i<this->getLength()) {
	std::forward_list<Mem> l = this->getMems(i);
	std::cout << i << std::endl;
	for(std::forward_list<Mem>::iterator it1=l.begin(); it1 != l.end(); ++it1) {
	    Mem m = (*it1);
	    std::cout << m.toStr() << std::endl;
	}
	i++;
    }
}
