#include <iostream>
#include <bitset>
#include "/home/dylan/Documents/code/ntHash-AVX512-rs_avx/ntHashIterator.hpp"
//#include "external/ntHash-AVX512/ntHashIterator.hpp"
#include <string>
#include <vector>

using namespace std;

// --------------------------------------------------

template<int BITS> class uint_kmer {
   using int_type =
      typename std::conditional< BITS <= 8,  uint8_t,
      typename std::conditional< BITS <= 16, uint16_t,
      typename std::conditional< BITS <= 32, uint32_t, 
      typename std::conditional< BITS <= 64, uint64_t,
      __uint128_t >::type >::type >::type >::type;
   public:
   int_type i;

   uint_kmer(int a) : i(a){}

   int_type get() {return i;}
};

// --------------------------------------------------

class Kmer {
  public :
    std::vector<uint64_t> arr;
    ssize_t k;
    Kmer(ssize_t k_) : arr(1+((k_-1)/32)), k{k_} {}

    /*
    friend ostream& operator<<(ostream& out, const Kmer& kmer) {
        ssize_t i = 0;
        ssize_t k_ = kmer.k;

        while(k_>32){
            out << bitset<64>(kmer.arr.at(i)) << " ";
            i++;
            k_-32;
        }

        out << bitset<k_>(kmer.arr.at(i));

        return out;
    }
    */
};

// --------------------------------------------------

int main(){
    cout << "hello, i do nothing for now\n";

    cout << "----------------------------------------" << endl;

    /* test sequence */
	std::string seq = "GAGTGTCAAACATTCAGACAACAGCAGGGGTGCTCTGGAATCCTATGTGAGGAACAAACATTCAGGCCACAGTAG";
	
	/* k is the k-mer length */
	unsigned k = 70;

    string kmer = seq.substr(0, k);
    uint64_t hVal, fhVal=0, rhVal=0; // canonical, forward, and reverse-strand hash values
    hVal = NTC64(kmer.c_str(), k, fhVal, rhVal); // initial hash value
    //...
    for (size_t i = 0; i < seq.length() - k; i++) {
        hVal = NTC64(seq[i], seq[i+k], k, fhVal, rhVal); // consecutive hash values
        cout << hVal << endl;
        cout << fhVal << endl;
        cout << rhVal << endl;
        cout << "----------" << endl;
        // hVal is the smallest of the two between fhVal and rhVal !
    }

    Kmer k1(20);
    Kmer k2(64);
    Kmer k3(65);
    Kmer k4(250);
    cout << k1.arr.size() <<" "<< k2.arr.size() <<" "<< k3.arr.size() <<" "<< k4.arr.size() << endl;

    return 0;
}