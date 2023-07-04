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

void splitIntoFile(std::string outfile, std::size_t id, const std::size_t k, 
                 const std::size_t m, const std::size_t fifo_size, Kmer** fifos,
                 sem_t* emptys, sem_t* fulls);

#endif