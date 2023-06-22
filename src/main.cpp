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
    bool isRevComp;

    Kmer(size_t k_, bool isRevComp_) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, isRevComp{isRevComp_}{}

    Kmer(size_t k_, bool isRevComp_, string seq) : arr(1+((k_-1)/32), 0), k{k_}, len{0}, isRevComp{isRevComp_}{
        // NOTE : try to make this function faster perhaps
        for (size_t i=0; i<seq.length(); i++){
            addNucl(seq[i]);
        }
    }

    Kmer(const Kmer& km): arr(km.arr), k{km.k}, len{km.len}, isRevComp{km.isRevComp} {}

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

    // Our k-sized window as we read the text, a.k.a. our "current" k-mer
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
    uint64_t fhpos; // position of the forward minimizer in array
    uint64_t fhmin; // forward minimizer hash value
    uint64_t rhpos; // position of the reverse-strand minimizer in array
    uint64_t rhmin; // reverse-strand minimizer hash value
    uint64_t hpos; // position of the canonical minimizer in array
    uint64_t hmin; // canonical minimizer hash value

    // Loop counter (used for our circular arrays fhvalues and rhvalues)
    uint16_t counter;

    // Super-k-mer FIFOs (buckets)
    queue<Kmer> fifos[q];
        // NOTE : maybe make it queues of Kmer* instead ?

    // Last read char
    char c;

    // -------------------- Program --------------------

    do {
        // To start, we need to read all the k first characters 
        for (size_t i=0; i<k; i++){
            kmer_cur[i] = fr.next_char();
        }

        // Let's compute the first hash values using ntHash
        for (size_t i=0; i<k-m+1; i++){
            NTC64(kmer_cur[i], kmer_cur[i+m], m, fhVal, rhVal);
            fhvalues[i] = fhVal;
            rhvalues[i] = rhVal;
            cout << fhVal <<" "<< rhVal << endl; // TEST
        }

        // Let's find the first minimizer 
        fhpos = min_element(fhvalues, fhvalues+k-m+1) - fhvalues;
        rhpos = min_element(rhvalues, rhvalues+k-m+1) - rhvalues;
        fhmin = fhvalues[fhpos];
        rhmin = rhvalues[rhpos];

        cout << endl << fhmin <<" "<< rhmin << endl; // TEST
        
        if (fhmin < rhmin){ 
            hmin = fhmin;
            hpos = fhpos;
            isRevComp = false;
        } else { // note : if they are equal, we arbitrarily pick the forward one
            hmin = rhmin;
            hpos = rhpos;
            isRevComp = true;
        } 

        // We create the first super-k-mer with our current k-mer
        Kmer* sk = new Kmer(2*k-m, isRevComp, kmer_cur);

        cout << hpos << " " << hmin << endl; // TEST

        c = fr.next_char();
        counter = 0;

        for(int ii=0; ii<10; ii++){ // TEST
        // while(c != '\0'){ // \0 should be returned at the end of a sequence (not of file)

            // Get the next k-mer (rotate the string once leftwise, and replace last character)
            rotate(kmer_cur.begin(), kmer_cur.begin()+1, kmer_cur.end());
            kmer_cur[k-1] = c;

            // Get the hash value of the new m-mer (rightmost)
            NTC64(kmer_cur[k-m-1], kmer_cur[k-1], m, fhVal, rhVal);

            // Place the new hash values in their respective arrays
            // Thanks to 'counter', we can use the array as a circular array
            fhvalues[counter] = fhVal;
            rhvalues[counter] = rhVal;

            /*
             * Now that we read a new character and therefore advanced the k-sized window,
             * we have many things to check, to see if it's the end of the current super-k-mer or not.
             * The possibilites are :
             * 
             * (1) the previous minimizer fell out of the k-sized window (if was the leftmost minimizer)
             * (2) the new m-mer (rightmost) is a better minimizer than the previous one
             * (3) the new m-mer (rightmost) is not better and we keep the previous one
             */

            // Let's check if the prev minimizer is still part of the new k-mer.
            if (hpos>0){

                // TODO

            } // If the minimizer fell out of the window, we end the super-k-mer.
            else {
                // We need to add the super-k-mer to its correct queue.
                fifos[hmin%q].push(*sk);

                // Now we need to create a new super-k-mer. First we need its minimizer.
                fhpos = min_element(fhvalues, fhvalues+k-m+1) - fhvalues;
                rhpos = min_element(rhvalues, rhvalues+k-m+1) - rhvalues;
                fhmin = fhvalues[fhpos];
                rhmin = rhvalues[rhpos];

                if (fhmin < rhmin){ 
                    hmin = fhmin;
                    hpos = fhpos;
                    isRevComp = false;
                } else { // note : if they are equal, we arbitrarily pick the forward one
                    hmin = rhmin;
                    hpos = rhpos;
                    isRevComp = true;
                } 

                // We can now create the new super-k-mer.
                sk = new Kmer(2*k-m, isRevComp, kmer_cur);
            }

            c = fr.next_char();
            counter = (counter+1)%(k-m+1);
        }
    
    } while (false); // TEST ; for now we don't want to loop over the whole file
    //} while (c != EOF)


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

    Kmer k1(35, "tttcccgggaaagggagggagagagagagagagagggaaagggaaaggg");
    // 32 first: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    cout << bitset<64>(k1.arr[0]) << endl;
    cout << bitset<64>(k1.arr[1]) << endl;

    return 0;
}