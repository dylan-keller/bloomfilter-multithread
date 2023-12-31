#ifndef KMER_HPP
#define KMER_HPP

#include <iostream>
#include <vector>
#include <string>

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
 * - size_t k : (maximum) size of the k-mer.
 *              a super-k-mer with m-sized minimizer should have 
 *              this attribute set to 2*k-m.
 * - size_t len : how many letters are in this k-mer.
 *                useful because otherwise we can't know if 00 is an
 *                empty char or an A.
 * - std::size_t position : used to store the position of the k-mer relative to the
 *                          input file. a position of i means this k-mer's first
 *                          letter is the i-th letter of the input file.
 *                          This will be used for the querying.
 * - bool isRevComp : if this k-mer is a reverse complement or not.
 */
class Kmer {
  public :
    std::vector<uint64_t> arr;
    std::size_t k;
    std::size_t len;
    std::size_t position;
    bool isRevComp;

    //Kmer();
    Kmer(std::size_t k_, bool isRevComp_);
    Kmer(std::size_t k_, std::size_t pos, bool isRevComp_);
    Kmer(std::size_t k_, bool isRevComp_, std::string seq);
    Kmer(std::size_t k_, std::size_t pos, bool isRevComp_, std::string seq);
    Kmer(const Kmer& km);

    void addNucl(char c);

    std::string to_string() const;

    friend std::ostream& operator<<(std::ostream& out, const Kmer& kmer);
};

#endif