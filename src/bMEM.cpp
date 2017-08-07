#include "bMEM.hpp"

BackwardMEM::BackwardMEM(const std::string& text, const std::string& path) {
    std::string index_path = path + ".idx";
    sdsl::construct_im(cst, text, 1);
}

tCST::size_type BackwardMEM::backward_search(const tCST::csa_type& csa,
                                             const unsigned char& c,
                                             tCST::size_type& lb,
                                             tCST::size_type& rb) {
    size_type cc = csa.char2comp[c];
    if(cc==0){
        lb=rb=0;
        return 0;
    }
    size_type c_before_l = csa.wavelet_tree.rank(lb, c);
    size_type c_before_r = csa.wavelet_tree.rank(rb, c);
    lb = csa.C[cc]+c_before_l;
    rb = csa.C[cc]+c_before_r;
    return rb-lb;
}

void BackwardMEM::addMEM(const std::string& s2,
                         const size_type& c_,
                         const size_type& p2_,
                         const size_type& k) {
    size_type sa_k = cst.csa[k];
    if(p2_ == 0) {
        MEMs.push_front(Mem(int(sa_k)+1, int(p2_)+1, int(c_)));
    } else {
        size_type lf_k = cst.csa.lf[k];
        size_type cc2 = cst.csa.char2comp[s2[p2_-1]];
        if(!(cst.csa.C[cc2] <= lf_k and lf_k < cst.csa.C[cc2+1])) {
            MEMs.push_front(Mem(int(sa_k)+1, int(p2_)+1, int(c_)));
        }
    }
}

std::list<Mem> BackwardMEM::getMEMs(const std::string& read,
                                    const unsigned int& l) {
    MEMs.clear();
    assert(l>0);

    size_type p2 = read.size();
    size_type i = 0, j = cst.size();
    size_type c = 0;

    while(p2>0){
        path_type path;
        size_type lb = i, rb = j;
        backward_search(cst.csa, read[p2-1], lb, rb);
        while(lb != rb and p2 > 0) {
            ++c;
            if(c >= l) {
                path.push_back(path_item(c, lb, rb, p2-1));
            }
            i = lb; j = rb;
            --p2;
            backward_search(cst.csa, read[p2-1], lb, rb);
        }
        for(path_type::iterator it=path.begin(); it!=path.end(); ++it) {
            size_type c_ = it->c;
            size_type lb_ = it->lb, rb_ = it->rb;
            lb = lb_; rb = lb_;
            size_type p2_ = it->p;
            while(c_ >= l) {
                for(size_type k = lb_; k < lb; ++k) {
                    addMEM(read, c_, p2_, k);
                }
                for(size_type k = rb; k < rb_; ++k) {
                    addMEM(read, c_, p2_, k);
                }
                lb = lb_; rb = rb_;
                node_type p = cst.parent( cst.node(lb_, rb_-1) );
                c_  = cst.depth(p);
                lb_ = cst.lb(p);
                rb_ = cst.rb(p)+1;
            }
        }
        if(0==c) {
            --p2;
        } else {
            node_type p = cst.parent( cst.node(i, j-1) );
            c = cst.depth(p);
            i = cst.lb(p);
            j = cst.rb(p)+1;
        }
    }
    return MEMs;
}

void BackwardMEM::printMEMs() {
    for(Mem m : MEMs) {
        std::cout << m.toStr() << " ";
    }
    std::cout << std::endl;
}
