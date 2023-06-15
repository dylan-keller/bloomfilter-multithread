#ifndef KMERHANDLER_HPP
#define KMERHANDLER_HPP

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

class KmerHandler{
  public:
    std::size_t k;
    uint64_t current_kmer;
    uint64_t mask;

    KmerHandler(std::size_t k_);
    ~KmerHandler();

    std::uint64_t next_Kmer(char nucleotide);
    std::string uint_To_Kmer(uint64_t val, std::size_t k);
};

#endif