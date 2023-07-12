#ifndef SKMERSPLITTER_HPP
#define SKMERSPLITTER_HPP

#include <fstream>
#include <iostream>
#include <semaphore.h>
#include <string>
#include <sys/stat.h>

#include <chrono>
#include <thread>

#include "Kmer.hpp"
#include "/home/dylan/Documents/code/bloomfilter-multithread/external/bitmagic/src/bm.h"

uint32_t xorshift32(const std::string& str);

void splitIntoFile(std::string outfile, std::size_t id, const std::size_t k, 
                 const std::size_t fifo_size, Kmer** fifo,
                 bool split_skmer_into_kmers, sem_t* empty, sem_t* full);

void splitIntoBF(std::size_t id, const std::size_t k, const std::size_t fifo_size, 
                 const std::size_t bf_size, Kmer** fifo, bm::bvector<>* bf, sem_t* empty, sem_t* full);

void splitQueryBF(std::size_t id, const std::size_t k, const std::size_t fifo_size,
                const std::size_t bf_size, Kmer** fifo, bm::bvector<>* bf, sem_t* empty, sem_t* full);

#endif