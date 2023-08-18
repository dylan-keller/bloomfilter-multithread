#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <queue>
#include <semaphore.h>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "../external/ntHash-AVX512/ntHashIterator.hpp"
#include "../external/bitmagic/src/bm.h"
#include "../external/bitmagic/src/bmserial.h"
#include "FastaReader.hpp"
#include "Kmer.hpp"
#include "KmerAnswer.hpp"
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

void invalid_message_exit() {
  std::cout << "Expected format is : \n\n"
            << "./main data_in query_in result_out k m t\n\n"
            << "With - data_in : path to FASTA file to be queried\n"
            << "     - query_in : path to FASTA file containing the query\n"
            << "     - result_out : path to file where the output will be saved\n"
            << "     - k : k-mer size\n"
            << "     - m : minimizer size, inferior to k\n"
            << "     - t : number of threads\n";
  exit(1);
}

// ------------------------------------------------------------

int main(int argc, char *argv[]){

    if (argc == 7) {
        std::cout << "FAILED : Invalid number of arguments. ";
        invalid_message_exit();
    }

    const uint16_t k = (uint16_t)atoi(argv[4]);
    const uint16_t m = (uint16_t)atoi(argv[5]);
    //const uint16_t q = (uint16_t)atoi(argv[6]);

    // cout << k << endl;

    // invalid_message_exit();

    // ------------------------------------------------------------

    std::cout << "----------------------------------------" << endl;

    // const size_t k = 30;
    // const size_t m = 15;
    const size_t q = 3;
    const size_t fifo_size = 100;
    const size_t bf_size = 65536;
    const size_t nb_query_buffers = 4;
    const size_t size_query_buffers = 512;
    atomic<size_t> finished_thread_counter(0);

    std::size_t fifo_out_counter[q] = {0};
    size_t buffer_counter[nb_query_buffers] = {0};

    Kmer* fifos_in[q*fifo_size];
    KmerAnswer* fifos_out[q*fifo_size];

    bm::bvector<> bloom_filters[q];
    bm::bvector<> query_answers[nb_query_buffers];

    sem_t emptys_in[q];
    sem_t fulls_in[q];
    sem_t emptys_out[q];
    sem_t fulls_out[q];

    for(size_t i=0; i<q; i++){
        bloom_filters[i].resize(bf_size);
        bloom_filters[i].set_range(0, bf_size, false);
        sem_init(&(emptys_in[i]), 0, fifo_size);
        sem_init(&(fulls_in[i]), 0, 0);
    }
    for(size_t i=0; i<nb_query_buffers; i++){
        query_answers[i].resize(size_query_buffers);
        query_answers[i].set_range(0, size_query_buffers, false);
    }

    thread extractor_thread(extractSkmers, "../inputs/sars-cov-2.fasta",
            k, m, q, fifo_size, fifos_in, emptys_in, fulls_in, nullptr);

    thread splitter_threads[q];

    /*
    for(size_t i=0; i<q; i++){
        splitter_threads[i] = thread(splitIntoFile, "../testing/", i, k,
        //                              fifo_size, fifos_in, true, &emptys_in[i], &fulls_in[i]);
                                        fifo_size, fifos_in, false, &emptys_in[i], &fulls_in[i]);
    }
    */
    
    for(size_t i=0; i<q; i++){
        splitter_threads[i] = thread(splitIntoBF, i, k, fifo_size, bf_size, fifos_in, &bloom_filters[i], &emptys_in[i], &fulls_in[i]);
    }

    extractor_thread.join();

    for(size_t i=0; i<q; i++){
        splitter_threads[i].join();
    }

    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys_in[i]));
        sem_destroy(&(fulls_in[i]));
    }

    std::cout << "---------- construction  done ----------" << endl;
    std::cout << "------------- query  start -------------" << endl;

    for(size_t i=0; i<q; i++){
        sem_init(&(emptys_in[i]), 0, fifo_size);
        sem_init(&(fulls_in[i]), 0, 0);
        sem_init(&(emptys_out[i]), 0, fifo_size);
        sem_init(&(fulls_out[i]), 0, 0);
    }

    thread extractor_thread2(extractSkmers, "../inputs/query.txt",
            k, m, q, fifo_size, fifos_in, emptys_in, fulls_in, &finished_thread_counter);

    thread splitter_threads2[q];
    
    for(size_t i=0; i<q; i++){
        splitter_threads2[i] = thread(splitQueryBF, i, k, fifo_size, bf_size,
                                      size_query_buffers, nb_query_buffers, fifos_in,
                                      fifos_out, &bloom_filters[i], &emptys_in[i], 
                                      &fulls_in[i], &emptys_out[i], &fulls_out[i],
                                      &finished_thread_counter);
    }

    // --------------------------------------------------

    KmerAnswer* skanswer;
    int semvalue = 0;
    size_t blocked_buffer = nb_query_buffers-1; // last bitvector is blocked
    size_t next_blocked_buffer = 0;
    bool fifos_are_empty = false;

    std::ofstream outputfile;
    outputfile.open("../outputs/answers.bin", std::ios::binary | std::ios::app);
  
    while((finished_thread_counter != q+1) || (!fifos_are_empty)){ // while all other threads are not finished :
        // std::this_thread::sleep_for(std::chrono::milliseconds(50)); // VALGRIND doesn't work if there is no sleep for some reason 
        for(size_t i=0; i<q; i++){
            sem_getvalue(&fulls_out[i], &semvalue); // check how many elements are in thread #i 's fifo
            while (semvalue > 0) { // while there are still elements to read,
                // read oldest element
                skanswer = fifos_out[i*fifo_size + fifo_out_counter[i]];

                // if it is a super-k-mer that would be written in the blocked bitvector,
                // we ignore it for now and stop reading this thread
                if (skanswer->outbv_nb == blocked_buffer) break;

                // else, we write it in our output bitvector

                // skanswer->fillOutputVector(&query_answers[skanswer->outbv_nb]); // TODO : overlap
                skanswer->fillOutputVector(
                    &query_answers[skanswer->outbv_nb], 
                    &query_answers[(skanswer->outbv_nb+1)%nb_query_buffers], 
                    size_query_buffers, 
                    &buffer_counter[skanswer->outbv_nb],
                    &buffer_counter[((skanswer->outbv_nb+1)%nb_query_buffers)]
                );

                // since we finished reading a KmerAnswer, we :
                // update the counter,
                fifo_out_counter[i] = (fifo_out_counter[i]+1) % fifo_size;
                // delete the KmerAnswer,
                delete skanswer;
                // and update the semaphores to notify a new empty spot in the corresponding fifo
                sem_wait(&fulls_out[i]);
                sem_post(&emptys_out[i]);
                // we also update our semvalue to read the appropriate number of elements
                semvalue--;
            }
        }
        next_blocked_buffer = (blocked_buffer+1)%nb_query_buffers;

        std::cout << "<" << buffer_counter[next_blocked_buffer] << ">\n";
        while (buffer_counter[next_blocked_buffer] == size_query_buffers){
            // bm::serializer<bm::bvector<>>::serialize(
            //      query_answers[next_blocked_buffer], outputfile, size_query_buffers
            // );

            // TODO : write in output file 
            
            query_answers[next_blocked_buffer].set_range(0, size_query_buffers, false);
            buffer_counter[next_blocked_buffer] = 0;
            blocked_buffer = next_blocked_buffer;
            next_blocked_buffer = (blocked_buffer+1)%nb_query_buffers;
            std::cout << '/' << blocked_buffer << "/\n";
        }
        
        fifos_are_empty = true;
        for(size_t i=0; i<q; i++){
            sem_getvalue(&fulls_out[i], &semvalue);
            std::cout << semvalue;
            if (semvalue > 0){ 
                std::cout << "!a!\n";
                fifos_are_empty = false;
                break;
            }
            std::cout << "\n";
        }
    }

    // to finish threads cleanly (although they are already done when this code is reached)
    extractor_thread2.join();
    for(size_t i=0; i<q; i++) splitter_threads2[i].join();

    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys_in[i]));
        sem_destroy(&(fulls_in[i]));
        sem_destroy(&(emptys_out[i]));
        sem_destroy(&(fulls_out[i]));
    }

    std::cout << "--------------- all done ---------------\n";

    // visualizeBitVector(bloom_filters[0]);
    // visualizeBitVector(query_answers[0]);
    // cout << "counter 0: " << counter[0] << endl;

    visualizeBitVector(query_answers[0]);
    visualizeBitVector(query_answers[1]);
    visualizeBitVector(query_answers[2]);
    visualizeBitVector(query_answers[3]);
    std::cout << buffer_counter[0] << " ----- " << buffer_counter[1] << " ----- " ;
    std::cout << buffer_counter[2] << " ----- " << buffer_counter[3] << endl;
    

    outputfile.close();

    return 0;
}