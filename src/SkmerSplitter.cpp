#include "../include/SkmerSplitter.hpp"
#include "SkmerSplitter.hpp"

QueryParameters::QueryParameters(
        std::size_t id, const std::size_t k, const std::size_t fifo_size, 
        const std::size_t bf_size, const std::size_t outbv_size, 
        const std::size_t outbv_nb, 
        Kmer** fifo, bm::bvector<>* bf, bm::bvector<>* outbv, 
        std::atomic<std::size_t>* counters, 
        sem_t* empty, sem_t* full) : 
        id{id}, k{k}, fifo_size{fifo_size}, bf_size{bf_size}, 
        outbv_size{outbv_size}, outbv_nb{outbv_nb},
        fifo{fifo}, bf{bf}, outbv{outbv}, counters{counters}, empty{empty},
        full{full} {}

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
                  const std::size_t outbv_size, const std::size_t outbv_nb, std::atomic<std::size_t>* counter, 
                  Kmer** fifo, bm::bvector<>* bf, bm::bvector<>* outbv, std::mutex* query_mutex, sem_t* empty, sem_t* full){
    std::size_t fifo_counter = 0;
    std::size_t counter_c = 0;
    std::size_t last_counter = 0;
    std::size_t kmer_pos = 0;
    std::size_t expected = 0;
    std::size_t desired;

    Kmer* sk;

    std::cout << "[query thread " << id << " start]" << std::endl;

    while(true){
        sem_wait(full);

        sk = fifo[id*fifo_size + fifo_counter];
        fifo_counter = (fifo_counter+1)%fifo_size;

        if((*sk).len == 0){ // sent if fasta file is finished
            counter[outbv_nb]++;
            std::cout << "[query thread " << id << " over]\n";
            delete sk;
            return;
        
        } else {
            kmer_pos = (*sk).position;
            std::cout << '{' << kmer_pos << '}' << (*outbv).test(2);
			std::string skstr = (*sk).to_string();

			for(std::size_t i=0; i<(*sk).len-k+1; i++){
                counter_c = ((kmer_pos+i)/outbv_size) % outbv_nb;

                while(counter_c != last_counter){
                    last_counter = (last_counter+1) % outbv_nb;
                    // check the lock number last_counter
                    {
                        std::cout << "<thread " << id << " will lock at " << last_counter << ">\n";
                        std::lock_guard<std::mutex> lock(query_mutex[last_counter]);
                        std::cout << "<thread " << id << " continues>\n";
                    }
                }
            
                uint32_t hash = (xorshift32(skstr.substr(i,k))) % bf_size;
                
                if ((*bf).test(hash)){
                    // (*outbv).set((kmer_pos+i) % outbv_size);
                    (*outbv)[(kmer_pos+i) % outbv_size] = true;
                    // if ((id==2)&&(kmer_pos<50)) std::cout << '<' << (kmer_pos+i) << '|' << (*outbv).test(kmer_pos+i) << '>';
                } 

                do {
                    desired = expected+1;
                } while(!counter[counter_c].compare_exchange_weak(expected, desired, std::memory_order_relaxed));

                // kmer_pos++;
			}
			delete sk;
        }
        sem_post(empty);
    }
}

/*
void splitQueryBF(QueryParameters qp){
    std::size_t fifo_counter = 0;
    std::size_t counter_c = 0;
    std::size_t kmer_pos = 0;
    std::size_t expected = 0;
    std::size_t desired;

    Kmer* sk;

    std::cout << "[query thread " << qp.id << " start]" << std::endl;

    while(true){
        sem_wait(qp.full);

        sk = qp.fifo[qp.id*qp.fifo_size + fifo_counter];
        fifo_counter = (fifo_counter+1)%qp.fifo_size;

        if((*sk).len == 0){ // sent if fasta file is finished
            std::cout << "[query thread " << qp.id << " over]\n";
            delete sk;
            return;
        
        } else {
            kmer_pos = (*sk).position;
            counter_c = (kmer_pos/qp.fifo_size) % qp.outbv_nb;

			std::string skstr = (*sk).to_string();
			for(std::size_t i=0; i< (*sk).len-qp.k+1; i++){
                
				(*qp.bf).set( (xorshift32(skstr.substr(i,qp.k)))%qp.bf_size );

                do {
                    desired = expected+1;
                } while(!qp.counters[counter_c].compare_exchange_weak(expected, desired, std::memory_order_relaxed));

                kmer_pos++;
			}
			delete sk;
        }
        sem_post(qp.empty);
    }
}
*/