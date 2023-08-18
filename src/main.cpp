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

    std::cout << "----------------------------------------" << endl;

    // -------------------- CREATING VARIABLES --------------------

    const size_t q = 3;                         // number of splitter threads
    const size_t fifo_size = 100;               // size of each thread's fifo 
    const size_t bf_size = 65536;               // size of each bloom filter
    const size_t nb_query_buffers = 4;          // number of buffers for the output
    const size_t size_query_buffers = 512;      // size of buffers for the output
    atomic<size_t> finished_thread_counter(0);  // will count how many query threads are done

    std::size_t fifo_out_counter[q] = {0};         // denotes which slot of the fifo to check for each thread
    size_t buffer_counter[nb_query_buffers] = {0}; // saves how many bits have been given to each output buffer

    Kmer* fifos_in[q*fifo_size];        // fifo for the super-k-mers given to the splitter threads
    KmerAnswer* fifos_out[q*fifo_size]; // fifo for the answers given by the splitter threads

    bm::bvector<> bloom_filters[q];                // bloom filters bit vectors
    bm::bvector<> query_answers[nb_query_buffers]; // buffers for the output

    sem_t emptys_in[q];  // semaphore denoting how many empty spots are in each input fifos
    sem_t fulls_in[q];   // semaphore denoting how many full spots are in each input fifos
    sem_t emptys_out[q]; // semaphore denoting how many empty spots are in each output fifos
    sem_t fulls_out[q];  // semaphore denoting how many full spots are in each output fifos

    // -------------------- PREPARING VARIABLES --------------------

    for(size_t i=0; i<q; i++){
        // setting up the bloom filters at the correct size and they start full of 0
        bloom_filters[i].resize(bf_size);
        bloom_filters[i].set_range(0, bf_size, false);
        // setting up the semaphores at the correct values
        sem_init(&(emptys_in[i]), 0, fifo_size);
        sem_init(&(fulls_in[i]), 0, 0);
    }

    // preparing the output bitvectors for the querying
    for(size_t i=0; i<nb_query_buffers; i++){
        query_answers[i].resize(size_query_buffers);
        query_answers[i].set_range(0, size_query_buffers, false);
    }

    // -------------------- STARTING FIRST SET OF THREADS --------------------
    // these threads will read the first input and use it to fill the bloom filters

    thread extractor_thread(extractSkmers, "inputs/sars-cov-2.fasta",
            k, m, q, fifo_size, fifos_in, emptys_in, fulls_in, nullptr);

    thread splitter_threads[q];
    
    for(size_t i=0; i<q; i++){
        splitter_threads[i] = thread(splitIntoBF, i, k, fifo_size, bf_size, fifos_in, &bloom_filters[i], &emptys_in[i], &fulls_in[i]);
    }
    
    // we wait for the extractor thread to finish first, then all splitters
    extractor_thread.join();
    for(size_t i=0; i<q; i++){
        splitter_threads[i].join();
    }

    // we destroy these semaphores since we are done with the first set of fifos
    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys_in[i]));
        sem_destroy(&(fulls_in[i]));
    }

    std::cout << "--------- bloom filters filled ---------" << endl;
    std::cout << "------------- query  start -------------" << endl;

    // -------------------- STARTING SECOND SET OF THREADS --------------------
    // these threads will query the bloom filters while the main saves their answers

    // setting up the semaphores at the correct values
    for(size_t i=0; i<q; i++){
        sem_init(&(emptys_in[i]), 0, fifo_size);
        sem_init(&(fulls_in[i]), 0, 0);
        sem_init(&(emptys_out[i]), 0, fifo_size);
        sem_init(&(fulls_out[i]), 0, 0);
    }

    thread extractor_thread2(extractSkmers, "inputs/query.txt",
            k, m, q, fifo_size, fifos_in, emptys_in, fulls_in, &finished_thread_counter);

    thread splitter_threads2[q];
    
    for(size_t i=0; i<q; i++){
        splitter_threads2[i] = thread(splitQueryBF, i, k, fifo_size, bf_size,
                                      size_query_buffers, nb_query_buffers, fifos_in,
                                      fifos_out, &bloom_filters[i], &emptys_in[i], 
                                      &fulls_in[i], &emptys_out[i], &fulls_out[i],
                                      &finished_thread_counter);
    }
    
    // -------------------- RECEIVING THREADS' ANSWERS --------------------
    // the set of threads are running, the main will now receive their answers and write them for outputting

    // ---------- setting some variables

    KmerAnswer* skanswer;                       // saves the last answer received 
    int semvalue = 0;                           // saves the value of the last semaphore checked

    // at all times there is always 1 blocked bitvector to prevent overwriting data that hasn't been saved
    size_t blocked_buffer = nb_query_buffers-1; // saves the number of the blocked bitvector
    size_t next_blocked_buffer = 0;             // saves the number of the next blocked bitvector

    bool fifos_are_empty = false;               // is only set to 'true' when all fifos are empty

    // ---------- setting an output stream

    std::ofstream outputfile;
    outputfile.open("../outputs/answers.bin", std::ios::binary | std::ios::app);

    // ---------- starting loop, will run until all threads are finished and all fifos are empty
  
    while((finished_thread_counter != q+1) || (!fifos_are_empty)){
        // std::this_thread::sleep_for(std::chrono::milliseconds(50)); // used for debugging if necessary

        // for each output fifo :
        for(size_t i=0; i<q; i++){
            sem_getvalue(&fulls_out[i], &semvalue); // check how many elements are in splitter thread #i 's fifo
            while (semvalue > 0) { // while there are still elements to read,
                // read oldest element
                skanswer = fifos_out[i*fifo_size + fifo_out_counter[i]];

                // if it is a super-k-mer that would be written in the blocked bitvector,
                // we ignore it for now and stop reading this thread's fifo
                if (skanswer->outbv_nb == blocked_buffer) break;

                // else, we write it in our output bitvector

                // this function will fill the correct output bitvector, and take care of any overlaps
                skanswer->fillOutputVector(
                    &query_answers[skanswer->outbv_nb], 
                    &query_answers[(skanswer->outbv_nb+1)%nb_query_buffers], 
                    size_query_buffers, 
                    &buffer_counter[skanswer->outbv_nb],
                    &buffer_counter[((skanswer->outbv_nb+1)%nb_query_buffers)]
                );

                // since we finished reading a KmerAnswer, we :
                // update the counter that tracks which slot of the fifo to read,
                fifo_out_counter[i] = (fifo_out_counter[i]+1) % fifo_size;
                // delete the KmerAnswer,
                delete skanswer;
                // and update the semaphores to notify there is a new empty spot in the corresponding fifo
                sem_wait(&fulls_out[i]);
                sem_post(&emptys_out[i]);
                // we also update our semvalue to read the appropriate number of elements
                semvalue--;
            }
        }

        // when we checked all fifos once, we see if we need to update which bitvector is blocked or not

        // (this variable is used to avoid doing this calculation multiple times 
        // and making the code easier to understand)
        next_blocked_buffer = (blocked_buffer+1)%nb_query_buffers;

        // we check if the oldest non-blocked bitvector is full
        // it's a 'while' and not an 'if' so that if it is indeed full, we also check the next oldest one
        while (buffer_counter[next_blocked_buffer] == size_query_buffers){
            // if it is full, we need to write it in memory, then clear the bitvector,
            // then unblock the previous one and block this one instead

            /*
            bm::serializer<bm::bvector<>>::serialize(
                 query_answers[next_blocked_buffer], outputfile, size_query_buffers
            );
            */

            // TODO : write in output file !!!
            
            // we clear the bitvector
            query_answers[next_blocked_buffer].set_range(0, size_query_buffers, false);
            // we set its counter back to zero
            buffer_counter[next_blocked_buffer] = 0;
            // we unblock the blocked buffer and block the next one instead
            blocked_buffer = next_blocked_buffer;
            // and we update the value of the 'next_blocked_buffer' (just for the while loop)
            next_blocked_buffer = (blocked_buffer+1)%nb_query_buffers;
        }
        
        // lastly we need to check if all fifos are empty or not
        fifos_are_empty = true;
        for(size_t i=0; i<q; i++){
            sem_getvalue(&fulls_out[i], &semvalue);
            if (semvalue > 0){ // semvalue is 0 iff the fifo is empty
                fifos_are_empty = false;
                break;
            }
        }
    }

    // to finish threads cleanly (although they are already done when this code is reached)
    extractor_thread2.join();
    for(size_t i=0; i<q; i++) splitter_threads2[i].join();

    // we destroy the semaphores since we are done with the fifos
    for(size_t i=0; i<q; i++){
        sem_destroy(&(emptys_in[i]));
        sem_destroy(&(fulls_in[i]));
        sem_destroy(&(emptys_out[i]));
        sem_destroy(&(fulls_out[i]));
    }

    std::cout << "--------------- all done ---------------\n";

    visualizeBitVector(query_answers[0]);
    visualizeBitVector(query_answers[1]);
    visualizeBitVector(query_answers[2]);
    visualizeBitVector(query_answers[3]);
    std::cout << buffer_counter[0] << " ----- " << buffer_counter[1] << " ----- " ;
    std::cout << buffer_counter[2] << " ----- " << buffer_counter[3] << endl;
    

    outputfile.close();

    return 0;
}