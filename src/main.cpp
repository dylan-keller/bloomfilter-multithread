#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <semaphore.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// #include "/mnt/c/Users/dylan/Documents/M1/stage/projet2/bloomfilter-multithread/external/ntHash-AVX512/ntHashIterator.hpp"
// #include "/mnt/c/Users/dylan/Documents/M1/stage/projet2/bloomfilter-multithread/external/bitmagic/src/bm.h"
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

    const size_t k = 30;
    const size_t m = 15;
    const size_t q = 3;
    const size_t fifo_size = 100;
    const size_t bf_size = 1024;
    const size_t nb_query_buffers = 8;
    const size_t size_query_buffers = 1024;
    size_t next_locked_buffer = 0;

    Kmer* fifos[q*fifo_size];

    bm::bvector<> bfs[q];

    atomic<size_t> counter[nb_query_buffers+1] = {}; // = {} initializes all values to zero
        // last counter used to count finished threads
    bm::bvector<> query_answers[nb_query_buffers];
    mutex query_mutex[nb_query_buffers];

    query_mutex[nb_query_buffers-1].lock(); // we lock the last buffer

    sem_t emptys[q];
    sem_t fulls[q];
    for(size_t i=0; i<q; i++){
        bfs[i].resize(bf_size);
        bfs[i].set_range(0, bf_size, false);
        query_answers[i].resize(size_query_buffers);
        query_answers[i].set_range(0, size_query_buffers, false);
        sem_init(&(emptys[i]), 0, fifo_size);
        sem_init(&(fulls[i]), 0, 0);
    }

    thread extractor_thread(extractSkmers, "/home/dylan/Documents/sequences/sars-cov-2.fasta",
            k, m, q, fifo_size, fifos, emptys, fulls, nullptr);

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

    cout << "---------- construction  done ----------" << endl;
    cout << "------------- query  start -------------" << endl;

    for(size_t i=0; i<q; i++){
        sem_init(&(emptys[i]), 0, fifo_size);
        sem_init(&(fulls[i]), 0, 0);
    }

    thread extractor_thread2(extractSkmers, "/home/dylan/Documents/sequences/query.txt",
            k, m, q, fifo_size, fifos, emptys, fulls, &counter[nb_query_buffers]);

    thread splitter_threads2[q];
    
    for(size_t i=0; i<q; i++){
        splitter_threads2[i] = thread(splitQueryBF, i, k, fifo_size, bf_size, 
                                    size_query_buffers, nb_query_buffers, counter, 
                                    fifos, &bfs[i], query_answers, 
                                    query_mutex, &emptys[i], &fulls[i]);
    }
  
    while(counter[nb_query_buffers] != q+1){
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        cout << "<new loop " << counter[nb_query_buffers] << ">\n";

        if (counter[next_locked_buffer] == size_query_buffers){
            query_mutex[next_locked_buffer].lock();
            counter[next_locked_buffer] = 0;
            query_mutex[(next_locked_buffer + nb_query_buffers-1)%nb_query_buffers].unlock(); //unlock the previous one

            next_locked_buffer = (next_locked_buffer+1)%nb_query_buffers;
        }
    }

    extractor_thread2.join();

    for(size_t i=0; i<q; i++){
        splitter_threads2[i].join();
    }

    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys[i]));
        sem_destroy(&(fulls[i]));
    }

    cout << "--------------- all done ---------------\n";

    //visualizeBitVector(bfs[0]);

    visualizeBitVector(query_answers[0]);

    std::this_thread::sleep_for(std::chrono::milliseconds(200));

    visualizeBitVector(query_answers[0]);

    cout << "counter 0: " << counter[0] << endl;

    std::cout << query_answers[0].test(1) << query_answers[0].test(2) << query_answers[0].test(3) << endl;

    return 0;
}