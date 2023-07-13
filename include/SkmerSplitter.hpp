#ifndef SKMERSPLITTER_HPP
#define SKMERSPLITTER_HPP

#include <atomic>
#include <fstream>
#include <iostream>
#include <semaphore.h>
#include <string>
#include <sys/stat.h>

#include <chrono>
#include <thread>

#include "Kmer.hpp"
#include "/home/dylan/Documents/code/bloomfilter-multithread/external/bitmagic/src/bm.h"

struct QueryParameters{
    std::size_t id;
    const std::size_t k;
    const std::size_t fifo_size;
    const std::size_t bf_size;
    const std::size_t outbv_size;
    const std::size_t outbv_nb; 
    Kmer** fifo;
    bm::bvector<>* bf;
    bm::bvector<>* outbv;
    std::atomic<std::size_t>* counters;
    sem_t* empty;
    sem_t* full;
    
    QueryParameters(
        std::size_t id, const std::size_t k, const std::size_t fifo_size, 
        const std::size_t bf_size, const std::size_t outbv_size, 
        const std::size_t outbv_nb, Kmer** fifo, bm::bvector<>* bf, 
        bm::bvector<>* outbv, std::atomic<std::size_t>* counters, 
        sem_t* empty, sem_t* full);
};

uint32_t xorshift32(const std::string& str);

void splitIntoFile(std::string outfile, std::size_t id, const std::size_t k, 
                 const std::size_t fifo_size, Kmer** fifo,
                 bool split_skmer_into_kmers, sem_t* empty, sem_t* full);

void splitIntoBF(std::size_t id, const std::size_t k, const std::size_t fifo_size,
                 const std::size_t bf_size, Kmer** fifo, bm::bvector<>* bf, sem_t* empty, sem_t* full);

void splitQueryBF(std::size_t id, const std::size_t k, const std::size_t fifo_size, const std::size_t bf_size, 
                  std::atomic<std::size_t>& counter, Kmer** fifo, bm::bvector<>* bf, sem_t* empty, sem_t* full);

void splitQueryBF(QueryParameters qp);

#endif