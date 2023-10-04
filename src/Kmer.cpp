#include "../include/Kmer.hpp"

// ------------------------------ constructors ------------------------------

Kmer::Kmer(std::size_t k_, bool isRevComp_) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, 
                                              position{0}, isRevComp{isRevComp_}{}

Kmer::Kmer(std::size_t k_, std::size_t pos, bool isRevComp_) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, 
                                                               position{pos}, isRevComp{isRevComp_}{}

Kmer::Kmer(std::size_t k_, bool isRevComp_, std::string seq) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, 
                                                               position{0}, isRevComp{isRevComp_}
{
    // NOTE : try to make this part faster perhaps
    for (std::size_t i=0; i<seq.length(); i++) addNucl(seq[i]);
}

Kmer::Kmer(std::size_t k_, std::size_t pos, bool isRevComp_, std::string seq) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, 
                                                                                position{pos}, isRevComp{isRevComp_}
{
    for (std::size_t i=0; i<seq.length(); i++) addNucl(seq[i]);
}

Kmer::Kmer(const Kmer& km): arr(km.arr), k{km.k}, len{km.len}, isRevComp{km.isRevComp} {}

// ------------------------------ functions ------------------------------

void Kmer::addNucl(char c){
    arr.at(len/32) |= ((c>>1)&(3UL)) << (2*(len%32));
    len++;
}

std::string Kmer::to_string() const{
    std::string res;
    char c;
    res.resize(len);

    if (!isRevComp){
        for(std::size_t i=0; i<len; i++){
            switch((arr.at(i/32) >> (2*(i%32))) & 3UL) {
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
        for(std::size_t i=len; i>0; i--){
            switch((arr.at((i-1)/32) >> (2*((i-1)%32))) & 3UL) {
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
            res[len-i] = c;
        }
    }
    return res;
}

std::ostream& operator<<(std::ostream& out, const Kmer& kmer) {
    out << kmer.to_string();
    return out;
}