#ifndef SKMERSPLITTER_HPP
#define SKMERSPLITTER_HPP

#include <atomic>
#include <fstream>
#include <iostream>
#include <mutex>
#include <semaphore.h>
#include <string>
#include <sys/stat.h>

#include <chrono>
#include <thread>

#include "Kmer.hpp"
#include "KmerAnswer.hpp"
#include "/mnt/c/Users/dylan/Documents/M1/stage/projet2/bloomfilter-multithread/external/bitmagic/src/bm.h"
// #include "/home/dylan/Documents/code/bloomfilter-multithread/external/bitmagic/src/bm.h"

uint32_t xorshift32(const std::string& str);

void splitIntoFile(std::string outfile, std::size_t id, const std::size_t k, 
                 const std::size_t fifo_size, Kmer** fifo,
                 bool split_skmer_into_kmers, sem_t* empty, sem_t* full);

void splitIntoBF(std::size_t id, const std::size_t k, const std::size_t fifo_size,
                 const std::size_t bf_size, Kmer** fifo, bm::bvector<>* bf, sem_t* empty, sem_t* full);

// I chose to keep outbv_size & outbv_nb as parameters, because that means each thread will do the modulo (%) operation on the k-mer position, 
// taking a little operation off the single thread filling the output (who is already the bottleneck).
void splitQueryBF(std::size_t id, const std::size_t k, const std::size_t fifo_size, const std::size_t bf_size,
                  const std::size_t outbv_size, const std::size_t outbv_nb, Kmer** fifo_in, KmerAnswer** fifo_out,
                  bm::bvector<>* bf, sem_t* empty_in, sem_t* full_in, sem_t* empty_out, sem_t* full_out,
                  std::atomic<std::size_t>* end_increment = nullptr);

#endif