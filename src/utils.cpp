#include "utils.hpp"

// TStr toTStr(const std::string& s) {
//     return TStr(s.c_str());
// }

int e_distance(const std::string& s1, const std::string& s2) {
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<int> col(len2+1), prevCol(len2+1);
    for(unsigned int i = 0; i < prevCol.size(); i++)
        prevCol[i] = i;
    for(unsigned int i = 0; i < len1; i++) {
        col[0] = i+1;
        for (unsigned int j = 0; j < len2; j++) {
            int m = (prevCol[1 + j] + 1)<(col[j] + 1)?(prevCol[1 + j] + 1):(col[j] + 1);
            col[j+1] = (m)<(prevCol[j] + (s1[i]==s2[j] ? 0 : 1))?(m):(prevCol[j] + (s1[i]==s2[j] ? 0 : 1));
        }
        col.swap(prevCol);
    }
    return prevCol[len2];
}

std::string reverse_and_complement(const std::string& s) {
    std::string rs = "";
    for(char c : s) {
        switch(c) {
        case 'A':
            rs += "T";
            break;
        case 'C':
            rs += "G";
            break;
        case 'G':
            rs += "C";
            break;
        case 'T':
            rs += "A";
            break;
        default:
            break;
        }
    }
    std::reverse(rs.begin(), rs.end());
    return rs;
}

bool compareMEMs(const Mem& m1, const Mem& m2) {
    if(m1.t <= m2.t)
        return false;
    else
        return true;
}
