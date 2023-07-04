#include <algorithm>
#include <chrono>
#include <iostream>
#include <queue>
#include <semaphore.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "/home/dylan/Documents/code/ntHash-AVX512-rs_avx/ntHashIterator.hpp"
#include "/home/dylan/Documents/code/bloomfilter-multithread/external/bitmagic/src/bm.h"
//#include "external/ntHash-AVX512/ntHashIterator.hpp"
#include "FastaReader.hpp"
#include "Kmer.hpp"
#include "SkmerExtractor.hpp"
#include "SkmerSplitter.hpp"

using namespace std;

// ------------------------------------------------------------

void visualizeBitVector(const bm::bvector<>& bv) {
    for (bm::id_t i = 0; i < bv.size(); ++i) {
        std::cout << (bv.test(i) ? "1" : "0");
        
        // Add line breaks for better readability
        if ((i + 1) % 80 == 0) {
            std::cout << std::endl;
        }
    }
    
    // Add a line break at the end if needed
    if (bv.size() % 80 != 0) {
        std::cout << std::endl;
    }
}

// ------------------------------------------------------------

void threadfun(queue<Kmer>* fifo, int i, sem_t* empty, sem_t* full){
    //this_thread::sleep_for(std::chrono::milliseconds(500*(i+1)));
    while(true){
        sem_wait(full);
        if(fifo->front().len == 0){ // sent if fasta file is finished
            cout << "[thread " << i << " over]\n";
            return;
        } else {
            cout << "(" << i << " " << fifo->front() << " " <<")\n";
            fifo->pop();
        }
        sem_post(empty);
    }
}

// ------------------------------------------------------------

int main(){
    cout << "hello, i do nothing for now\n";

    cout << "----------------------------------------" << endl;

    const size_t k = 12;
    const size_t m = 6;
    const size_t q = 10;
    const size_t fifo_size = 10;
    ssize_t id = -1;

    Kmer* fifos[q*fifo_size];

    sem_t emptys[q];
    sem_t fulls[q];
    for(size_t i=0; i<q; i++){
        sem_init(&(emptys[i]), 0, fifo_size);
        sem_init(&(fulls[i]), 0, 0);
    }

    // NOTE : these objects have the same size (might be worth noting for cache uses)
    // Kmer k1(10,false);
    // Kmer k2(1000,false);
    // cout << sizeof(k1) << " " << sizeof(k2) << endl;

    cout << "----------------------------------------" << endl;

    //extractSkmers("/home/dylan/Documents/sequences/sars-cov-2.fasta", 
    //              k, m, q, fifo_size, fifos, emptys, fulls);

    //splitIntoFile("../testing.txt", 0, k, m, fifo_size, fifos, &emptys[0], &fulls[0]);

    thread t1(extractSkmers, "/home/dylan/Documents/sequences/sars-cov-2.fasta",
            k, m, q, fifo_size, fifos, emptys, fulls);

    thread extractor_threads[q];

    for(size_t i=0; i<q; i++){
        extractor_threads[i] = thread(splitIntoFile, "../testing/testingX.txt", i, k, m,
                                      fifo_size, fifos, true, &emptys[i], &fulls[i]);
    }

    t1.join();

    for(size_t i=0; i<q; i++){
        extractor_threads[i].join();
    }

    cout << "----------------------------------------" << endl;

    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys[i]));
        sem_destroy(&(fulls[i]));
    }

    std::cout << "[------------]\n[- all done -]\n[------------]" << std::endl;

    return 0;
}