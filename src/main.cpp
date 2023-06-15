#include <iostream>
#include <bitset>
#include "/home/dylan/Documents/code/ntHash-AVX512-rs_avx/ntHashIterator.hpp"
//#include "external/ntHash-AVX512/ntHashIterator.hpp"
#include <string>

using namespace std;

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

int main(){
    cout << "hello, i do nothing for now\n";
    
    uint64_t a = 1ul;

    for (int k=1; k<33;){
        a = ((1ul << (2 * (k-1))) - 1);
        cout << k << " " << bitset<64>(a) << " " << a << endl;
        k = k*2;
    }

    cout << sizeof(uint8_t) << " " << sizeof(uint16_t) << " " << sizeof(uint32_t) << " " << sizeof(uint64_t) << endl;
    /*   
    uint_kmer<1> b1 = ((1ul << (20)) - 1); 
    uint_kmer<15> b2 = ((1ul << (20)) - 1);
    uint_kmer<16> b3 = ((1ul << (20)) - 1); 
    uint_kmer<17> b4 = ((1ul << (100)) -1);
    uint_kmer<99999> b5 = (unsigned long long)18446744073709551600+66;
    __uint128_t test128 = (unsigned long long)18446744073709551600+66;
    cout << sizeof(b1) <<" "<< sizeof(b2) <<" "<< sizeof(b3) <<" "<< sizeof(b4) <<" "<< sizeof(b5);
    cout << endl;
    cout << b1.get() <<" "<< b2.get() <<" "<< b4.get() <<" "<< (unsigned long long)b5.get() <<" "<< (unsigned long long)(test128);
    cout << endl;*/

    /* test sequence */
	std::string seq = "GAGTGTCAAACATTCAGACAACAGCAGGGGTGCTCTGGAATCCTATGTGAGGAACAAACATTCAGGCCACAGTAG";
	
	/* k is the k-mer length */
	unsigned k = 70;
	
	/* h is the number of hashes for each k-mer */
	unsigned h = 1;

	/* init ntHash state and compute hash values for first k-mer */
	ntHashIterator itr(seq, h, k);
	while (itr != itr.end()) {
		std::cout << (*itr)[0] << std::endl;
		++itr;
	}

	return 0;
}