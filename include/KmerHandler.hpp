#ifndef KMERHANDLER_HPP
#define KMERHANDLER_HPP

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

class KmerHandler{
  public:
    std::size_t k;
    __uint128_t current_kmer;
    __uint128_t mask;

    KmerHandler(std::size_t k_);
    ~KmerHandler();

    __uint128_t next_Kmer(char nucleotide);
    std::string uint_To_Kmer(__uint128_t val, std::size_t k);
};

#endif