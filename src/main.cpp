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

// ------------------------------------------------------------

/*
 * Kmer
 *
 * Class used to represent a k-mer or super-k-mer. Each nucleotide is represented
 * using 2 bits : 
 * A:00, C:01, T:10, G:11
 * 
 * Attributes :
 * - std::vector<uint64_t> arr : used to store the bits representing the k-mer. 
 *                               first letter is the 2 lowest bits of arr.at(0),
 *                               second letter is the 3rd and 4th bit of arr.at(0),
 *                               the 33rd letter is the 2 lowest bits of arr.at(1),
 *                               and so on.
 * - size_t k : size of the k-mer.
 *              a super-k-mer with m-sized minimizer should have 
 *              this attribute set to 2*k-m.
 * - size_t len : how many letters are in this k-mer.
 *                useful because otherwise we can't know if 00 is an
 *                empty char or an A.
 */
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

    // TODO : take care of reverse complement kmers (implement method or constructor)
    // or let the super-k-mer reader take care of it...

    friend ostream& operator<<(ostream& out, const Kmer& kmer) {
        string res;
        char c;
        res.resize(kmer.len);

        for(size_t i=0; i<kmer.len; i++){
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

        out << res;

        return out;
    }

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

    // Flag if the super-k-mer has ended and we need to create a new one
    bool new_skmer_flag = false;

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
        }

        // Let's find the first minimizer 
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

        // We create the first super-k-mer with our current k-mer
        Kmer* sk = new Kmer(2*k-m, isRevComp, kmer_cur);

        c = fr.next_char();
        counter = 0;

        // for(int ii=0; ii<60; ii++){ // TEST
        while(c != '\0'){ // \0 should be returned at the end of a sequence

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
            if (counter != hpos){ // If it is,
                // We need to check if the new m-mer (rightmost) is a better minimizer.
                if (fhvalues[counter] < hmin) {
                    hmin = fhvalues[counter];
                    hpos = counter;
                    isRevComp = false;
                    new_skmer_flag = true;
                } else if (rhvalues[counter] < hmin) {
                    hmin = rhvalues[counter];
                    hpos = counter;
                    isRevComp = true;
                    new_skmer_flag = true;
                } else {
                    sk->addNucl(c);
                }
            } 
            else { // If the minimizer fell out of the window, we end the super-k-mer.
                new_skmer_flag = true;

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
            }

            if (new_skmer_flag){
                // We need to add the super-k-mer to its correct queue.
                fifos[hmin%q].push(*sk);

                // We can now create the new super-k-mer, starting at the current k-mer.
                sk = new Kmer(2*k-m, isRevComp, kmer_cur);
            }

            c = fr.next_char();
            counter = (counter+1)%(k-m+1);
            new_skmer_flag = false;
        }
    
    } while (false); // TEST ; for now we don't want to loop over the whole file
    //} while (c != EOF)

    for (size_t i=0; i<q; i++){
        cout << fifos[i].size() << " " ;
        
        /*
        while(!fifos[i].empty()){
            cout << fifos[i].front() << endl;
            fifos[i].pop();
        }
        */
    }
    cout << endl;




}


// ------------------------------------------------------------

int main(){
    cout << "hello, i do nothing for now\n";

    cout << "----------------------------------------" << endl;

    run("/home/dylan/Documents/sequences/sars-cov-2.fasta", 129, 65, 5);

    cout << "----------------------------------------" << endl;

    return 0;
}