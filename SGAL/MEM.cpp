#include "MEM.hpp"

Mem::Mem() {
    this->t = 0;
    this->p = 0;
    this->l = 0;
}

Mem::Mem(int t, int p, int l) {
    this->t = t;
    this->p = p;
    this->l = l;
}

int Mem::getT() {
    return t;
}

int Mem::getP() {
    return p;
}

int Mem::getL() {
    return l;
}

std::string Mem::toStr() {
    return "(" + std::to_string(t) + "," +
	std::to_string(p) + "," +
	std::to_string(l) + ")";
}
