#include <iostream>
#include <bitset>
#include <string>
#include <vector>
#include <queue>
#include <stdexcept>

#include "/home/dylan/Documents/code/ntHash-AVX512-rs_avx/ntHashIterator.hpp"
//#include "external/ntHash-AVX512/ntHashIterator.hpp"
#include "FastaReader.hpp"
#include "KmerHandler.hpp"

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

// ------------------------------------------------------------

class Kmer {
  public :
    std::vector<uint64_t> arr;
    size_t k;
    size_t len;

    Kmer(size_t k_) : arr(1+((k_-1)/32), 0), k{k_}, len{0} {}

    Kmer(size_t k_, string seq) : arr(1+((k_-1)/32), 0), k{k_}, len{0} {
        // todo : make this function faster
        for (size_t i=0; i<seq.length(); i++){
            addNucl(seq[i]);
        }
    }

    /*
    friend ostream& operator<<(ostream& out, const Kmer& kmer) {
        size_t i = 0;
        size_t k_ = kmer.k;

        while(k_>32){
            out << bitset<64>(kmer.arr.at(i)) << " ";
            i++;
            k_-32;
        }

        out << bitset<k_>(kmer.arr.at(i));

        return out;
    }
    */

    void addNucl(char c){
        arr.at(len/32) |= ((c>>1)&(3UL)) << (2*(len%32));

        len++;
    }
};

// ------------------------------------------------------------

void run(string filename, const size_t k, const size_t m, const size_t q){

    // -------------------- Variables --------------------

    FastaReader fr(filename);
    KmerHandler kh(m);

    // Our k-sized window as we read the text
    string kmer_cur;
    kmer_cur.resize(k);

    // The forward (resp. reverse-strand) hash values of
    // every m-mer in the k-mer, for quick minimizer find
    uint64_t fhvalues[k-m+1];
    uint64_t rhvalues[k-m+1];

    // The forward (resp. reverse-strand) hash values
    uint64_t fhVal=0, rhVal=0; 

    // Tells us our kmer is reverse complement or not
    bool isRevComp;

    // The minimizer hash values :
    uint64_t* fhmin; // forward minimizer hash value
    uint64_t* rhmin; // reverse-strand minimizer hash value
    uint64_t hmin; // canonical minimizer hash value

    // Super-k-mer FIFOs (buckets)
    queue<Kmer> fifos[q];
        // note : maybe make it queues of Kmer* instead ?

    // -------------------- Program --------------------

    // At first, we need to read all the k first characters 
    for (size_t i=0; i<k; i++){
        kmer_cur[i] = fr.next_char();
    }

    // Let's compute the first hash values using ntHash
    for (size_t i=0; i<k-m+1; i++){
        NTC64(kmer_cur[i], kmer_cur[i+m], m, fhVal, rhVal);
        fhvalues[i] = fhVal;
        rhvalues[i] = rhVal;
    }

    // Let's find the first minimizer 
    fhmin = min_element(fhvalues, fhvalues+k-m+1);
    rhmin = min_element(rhvalues, rhvalues+k-m+1);
    
    if (*fhmin < *rhmin){ 
        hmin = *fhmin;
        isRevComp = false;
    } else { // note : if they are equal, we arbitrarily pick the forward one
        hmin = *rhmin;
        isRevComp = true;
    } 
    

    
}


// ------------------------------------------------------------

int main(){
    cout << "hello, i do nothing for now\n";

    cout << "----------------------------------------" << endl;

    run("/home/dylan/Documents/sequences/sars-cov-2.fasta", 20, 8, 10);
    cout << endl;

    cout << "----------------------------------------" << endl;

    string s = "aaabbbcccdddeee";
    rotate(s.begin(), s.begin() + 1, s.end());
    s[s.length()-1] = 'f';
    cout << s << endl;

    cout << "----------------------------------------" << endl;

    Kmer k1(60, "gagagagagagggagggagagagagagagagagggaaagggaaaggg");
    //           XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    //           gagagagagagggagggagagagagagagaga

    cout << bitset<64>(k1.arr[0]) << endl;
    cout << bitset<64>(k1.arr[1]) << endl;

    return 0;
}