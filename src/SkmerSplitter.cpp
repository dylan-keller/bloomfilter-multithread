#include "../include/SkmerSplitter.hpp"
#include "SkmerSplitter.hpp"

void splitSkmerIntoKmers();

void splitIntoFile(std::string outfile, std::size_t id, const std::size_t k,
                 const std::size_t m, const std::size_t fifo_size, Kmer **fifo,
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
		//std::this_thread::sleep_for(std::chrono::milliseconds(400*(id+1)));
        sem_wait(full);

		ssize_t truc = id*fifo_size + fifo_counter;

        sk = fifo[truc];
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