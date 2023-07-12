#include <algorithm>
#include <chrono>
#include <iomanip>
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

int main(){
    cout << "hello, i do nothing for now\n";

    cout << "----------------------------------------" << endl;

    const size_t k = 60;
    const size_t m = 30;
    const size_t q = 20;
    const size_t fifo_size = 100;
    const size_t bf_size = 1000;

    Kmer* fifos[q*fifo_size];

    bm::bvector<> bfs[q];

    sem_t emptys[q];
    sem_t fulls[q];
    for(size_t i=0; i<q; i++){
        bfs[i].resize(bf_size);
        bfs[i].set_range(0, bf_size, false);
        sem_init(&(emptys[i]), 0, fifo_size);
        sem_init(&(fulls[i]), 0, 0);
    }

    thread extractor_thread(extractSkmers, "/home/dylan/Documents/sequences/sars-cov-2.fasta",
            k, m, q, fifo_size, fifos, emptys, fulls);

    thread splitter_threads[q];

    /*
    for(size_t i=0; i<q; i++){
        splitter_threads[i] = thread(splitIntoFile, "../testing/", i, k,
        //                              fifo_size, fifos, true, &emptys[i], &fulls[i]);
                                        fifo_size, fifos, false, &emptys[i], &fulls[i]);
    }
    */
    
    for(size_t i=0; i<q; i++){
        splitter_threads[i] = thread(splitIntoBF, i, k, fifo_size, bf_size, fifos, &bfs[i], &emptys[i], &fulls[i]);
    }


    extractor_thread.join();

    for(size_t i=0; i<q; i++){
        splitter_threads[i].join();
    }

    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys[i]));
        sem_destroy(&(fulls[i]));
    }

    cout << "-------- construction  complete --------" << endl;

    /*

    for(size_t i=0; i<q; i++){
        sem_init(&(emptys[i]), 0, fifo_size);
        sem_init(&(fulls[i]), 0, 0);
    }

    thread extractor_thread2(extractSkmers, "/home/dylan/Documents/sequences/query.txt",
            k, m, q, fifo_size, fifos, emptys, fulls);

    thread splitter_threads2[q];

    for(size_t i=0; i<q; i++){
        splitter_threads2[i] = thread(splitIntoFile, "../testing/", i, k,
        //                              fifo_size, fifos, true, &emptys[i], &fulls[i]);
                                        fifo_size, fifos, false, &emptys[i], &fulls[i]);
    }

    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys[i]));
        sem_destroy(&(fulls[i]));
    }
    */

    cout << "--------------- all done ---------------\n";

    visualizeBitVector(bfs[0]);

    return 0;
}