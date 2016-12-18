#include "MEMsList.hpp"

MemsList::MemsList(const int& l) {
    length = l+1;
    std::vector<std::forward_list<Mem> > mems (length, std::forward_list<Mem> ());
    this->mems = mems;
}

void MemsList::addMem(const int& t, const int& p, const int& l) {
    mems[p].push_front(Mem(t,p,l));
}

std::forward_list<Mem> MemsList::getMems(const int& i) {
    return mems[i];
}

int MemsList::getLength() {
    return length;
}
