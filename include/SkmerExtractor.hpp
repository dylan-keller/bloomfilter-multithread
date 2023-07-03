#ifndef SKMEREXTRACTOR_HPP
#define SKMEREXTRACTOR_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <semaphore.h>
#include <thread>
#include <vector>

#include "FastaReader.hpp"
#include "Kmer.hpp"

#include "/home/dylan/Documents/code/ntHash-AVX512-rs_avx/ntHashIterator.hpp"

void extractSkmers(std::string filename, const std::size_t k, const std::size_t m,
                const std::size_t q, const std::size_t fifo_size, Kmer** fifos,
                sem_t* emptys, sem_t* fulls);

#endif