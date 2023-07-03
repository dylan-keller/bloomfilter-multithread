#include "../include/Kmer.hpp"

Kmer::Kmer(){}

Kmer::Kmer(std::size_t k_, bool isRevComp_) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, isRevComp{isRevComp_}{}

Kmer::Kmer(std::size_t k_, bool isRevComp_, std::string seq) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, isRevComp{isRevComp_}{
    // NOTE : try to make this function faster perhaps
    for (std::size_t i=0; i<seq.length(); i++){
        addNucl(seq[i]);
    }
}

Kmer::Kmer(const Kmer& km): arr(km.arr), k{km.k}, len{km.len}, isRevComp{km.isRevComp} {}


void Kmer::addNucl(char c){
    arr.at(len/32) |= ((c>>1)&(3UL)) << (2*(len%32));
    len++;
}

std::ostream& operator<<(std::ostream& out, const Kmer& kmer) {
    std::string res;
    char c;
    res.resize(kmer.len);

    if (!kmer.isRevComp){
        for(std::size_t i=0; i<kmer.len; i++){
            switch((kmer.arr.at(i/32) >> (2*(i%32))) & 3UL) {
                case 0:
                    c = 'A';
                    break;
                case 1:
                    c = 'C';
                    break;
                case 2:
                    c = 'T';
                    break;
                case 3:
                    c = 'G';
                    break;
            }
            res[i] = c;
        }
    } else {
        for(std::size_t i=kmer.len; i>0; i--){
            switch((kmer.arr.at((i-1)/32) >> (2*((i-1)%32))) & 3UL) {
                case 0:
                    c = 'T';
                    break;
                case 1:
                    c = 'G';
                    break;
                case 2:
                    c = 'A';
                    break;
                case 3:
                    c = 'C';
                    break;
            }
            res[kmer.len-i] = c;
        }
    }

    out << res;

    return out;
}
