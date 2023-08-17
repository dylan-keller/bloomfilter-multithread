#include "../include/SkmerSplitter.hpp"
#include "SkmerSplitter.hpp"

uint32_t xorshift32(const std::string& str) {
    uint32_t hash = 0;
    for (char c : str) {
        hash ^= c;
        hash ^= hash << 13;
        hash ^= hash >> 17;
        hash ^= hash << 5;
    }
    return hash;
}

void splitIntoFile(std::string outfile, std::size_t id, const std::size_t k,
                 const std::size_t fifo_size, Kmer **fifo,
                 bool split_skmer_into_kmers, sem_t *empty, sem_t *full) 
{
	if(outfile[outfile.length()-1] != '/')
		outfile += '/';
	outfile = outfile + std::to_string(id) + ".txt";

    std::ofstream outf(outfile);
    
    std::size_t fifo_counter = 0;
    Kmer* sk;

    std::cout << "[file thread " << id << " start]" << std::endl;

    while(true){
        sem_wait(full);

        sk = fifo[id*fifo_size + fifo_counter];
        fifo_counter = (fifo_counter+1)%fifo_size;

        if((*sk).len == 0){ // sent if fasta file is finished
            std::cout << "[file thread " << id << " over]\n";
            outf.close();
            delete sk;
            return;
        } else {
			std::string skstr = (*sk).to_string();
			/*
			for (auto i = skstr.cbegin(); i != skstr.cend() - k + 1; ++i) {
				std::string substring(&*i, k);
        		outf << substring << '\n';
			}*/
			if (split_skmer_into_kmers){
				for(std::size_t i=0; i< (*sk).len-k+1; i++){
					outf << skstr.substr(i,k) << '\n';
				}
			} else {
				outf << *sk << "\n";
			}
			delete sk;
        }
        sem_post(empty);
    }
}


void splitIntoBF(std::size_t id, const std::size_t k, const std::size_t fifo_size,
                 const std::size_t bf_size, Kmer** fifo, bm::bvector<>* bf, sem_t* empty, sem_t* full){
	std::size_t fifo_counter = 0;

    Kmer* sk;

    std::cout << "[fill thread " << id << " start]" << std::endl;

    while(true){
        sem_wait(full);

        sk = fifo[id*fifo_size + fifo_counter];
        fifo_counter = (fifo_counter+1)%fifo_size;

        if((*sk).len == 0){ // sent if fasta file is finished
            std::cout << "[fill thread " << id << " over]\n";
            delete sk;
            return;
        } else {
			std::string skstr = (*sk).to_string();
			for(std::size_t i=0; i< (*sk).len-k+1; i++){
				(*bf).set( (xorshift32(skstr.substr(i,k)))%bf_size );
			}
			delete sk;
        }
        sem_post(empty);
    }
}

void splitQueryBF(std::size_t id, const std::size_t k, const std::size_t fifo_size, const std::size_t bf_size,
                  const std::size_t outbv_size, const std::size_t outbv_nb, Kmer** fifo_in, KmerAnswer** fifo_out,
                  bm::bvector<>* bf, sem_t* empty_in, sem_t* full_in, sem_t* empty_out, sem_t* full_out,
                  std::atomic<std::size_t>* end_increment){
    std::size_t fifo_counter = 0;
    std::size_t fifo_arrayspot;

    Kmer* sk; // used to parse inputs
    KmerAnswer* skanswer; // used to create outputs

    std::cout << "[query thread " << id << " start]" << std::endl; // debugging

    while(true){
        sem_wait(full_in);   // we wait for a super-k-mer in the input fifo
        sem_wait(empty_out); // we wait for a free spot in the output fifo

        fifo_arrayspot = id*fifo_size + fifo_counter;
        sk = fifo_in[fifo_arrayspot];

        if((*sk).len == 0){ // sent iff fasta file is finished
            // counter[outbv_nb]++; needed to count how many threads have finished 
            std::cout << "[query thread " << id << " over]\n";
            if (end_increment != nullptr) end_increment[0]++; // TODO : check if atomic
            delete sk;
            return;
        
        } else { // we just recieved a super-k-mer we have to split
            // so we create the answer we will send back
            skanswer = new KmerAnswer(
                ((*sk).position)%outbv_size,             // starting position IN the output bitvector.
                (((*sk).position)/outbv_size) % outbv_nb  // which output bitvector it starts in.
            );
            skanswer->bv.resize((*sk).len-k+1);
            skanswer->bv.set_range(0, skanswer->bv.size(), false);
    
            // we get the super-k-mer string
			std::string skstr = (*sk).to_string();

            // for each k-mer of the super-k-mer :
			for(std::size_t i=0; i<(*sk).len-k+1; i++){
                // we hash the current k-mer
                uint32_t hash = (xorshift32(skstr.substr(i,k))) % bf_size;
                
                // we check if the current k-mer hash is present in this thread's bloom filter
                //if ((*bf).test(hash)) skanswer->bv.set(i);
                skanswer->bv.set(i, (*bf).test(hash));
			}
            // we finished parsing this super-k-mer, so we delete it
			delete sk;

            // and we send the answer in the fifo output, at the correct position
            fifo_out[fifo_arrayspot] = skanswer;

            // finally we increment the counter (for inputs and outputs)
            fifo_counter = (fifo_counter+1)%fifo_size;
        }
        sem_post(empty_in); // we notify that there is now a new empty spot in the input fifo
        sem_post(full_out); // we notify that there is now a filled spot in the output fifo
        
    }
}