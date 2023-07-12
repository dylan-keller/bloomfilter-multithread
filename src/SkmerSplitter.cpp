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

    std::cout << "[thread " << id << " start]" << std::endl;

    while(true){
        sem_wait(full);

        sk = fifo[id*fifo_size + fifo_counter];
        fifo_counter = (fifo_counter+1)%fifo_size;

        if((*sk).len == 0){ // sent if fasta file is finished
            std::cout << "[thread " << id << " over]\n";
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

    std::cout << "[thread " << id << " start]" << std::endl;

    while(true){
        sem_wait(full);

        sk = fifo[id*fifo_size + fifo_counter];
        fifo_counter = (fifo_counter+1)%fifo_size;

        if((*sk).len == 0){ // sent if fasta file is finished
            std::cout << "[thread " << id << " over]\n";
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

void splitQueryBF(std::size_t id, const std::size_t k, const std::size_t fifo_size,
                const std::size_t bf_size, Kmer** fifo, bm::bvector<>* bf, sem_t* empty, sem_t* full){
    std::size_t fifo_counter = 0;
    Kmer* sk;

    std::cout << "[thread " << id << " start]" << std::endl;

    while(true){
        sem_wait(full);

        sk = fifo[id*fifo_size + fifo_counter];
        fifo_counter = (fifo_counter+1)%fifo_size;

        if((*sk).len == 0){ // sent if fasta file is finished
            std::cout << "[thread " << id << " over]\n";
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