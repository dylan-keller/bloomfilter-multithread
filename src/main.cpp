#include <iostream>
#include <bitset>
#include <string>
#include <vector>

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

// ------------------------------------------------------------

void run(string filename, const size_t k, const size_t m){

    FastaReader fr(filename);
    KmerHandler kh(m);

    // Will represent our k-sized window as we read the text
    string kmer_cur;
    kmer_cur.resize(k);

    // At first, we need to read all the k first characters
    for (size_t i=0; i<k; i++){
        kmer_cur[i] = fr.next_char();
    }

    // Will be used to store the forward (resp. backward) hash values
    // of every m-mer in the k-mer, for quick minimizer find
    uint64_t fhvalues[k-m+1];
    uint64_t rhvalues[k-m+1];

    // Canonical, forward, and reverse-strand hash values
    uint64_t hVal, fhVal=0, rhVal=0; 

    // Let's compute the first hash values using ntHash
    for (size_t i=0; i<k-m+1; i++){
        NTC64(kmer_cur[i], kmer_cur[i+m], m, fhVal, rhVal);
        fhvalues[i] = fhVal;
        rhvalues[i] = rhVal;
        cout << fhvalues[i] << " " << rhvalues[i] << endl; // TEST
    }


    uint64_t* fhmin = min_element(fhvalues, fhvalues+k-m+1);
    uint64_t* rhmin = min_element(rhvalues, rhvalues+k-m+1);

    cout << endl << *fhmin << " " << *rhmin << endl;



    cout << endl << kmer_cur << endl; // TEST
    

}


// ------------------------------------------------------------

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
        cout << "---" << endl;
        // hVal is the smallest of the two between fhVal and rhVal !
    }

    Kmer k1(20);
    Kmer k2(64);
    Kmer k3(65);
    Kmer k4(250);
    cout << k1.arr.size() <<" "<< k2.arr.size() <<" "<< k3.arr.size() <<" "<< k4.arr.size() << endl;

    cout << "----------------------------------------" << endl;

    run("/home/dylan/Documents/sequences/sars-cov-2.fasta", 20, 8);

    cout << "----------------------------------------" << endl;

    string s = "aaabbbcccdddeee";
    rotate(s.begin(), s.begin() + 1, s.end());
    s[s.length()-1] = 'f';
    cout << s << endl;

    return 0;
}