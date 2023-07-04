#include "../include/SkmerSplitter.hpp"
#include "SkmerSplitter.hpp"

/*
class SkmerSaver
{
public:
	std::string outdir;
	std::unordered_map<uint64_t, std::stringstream> file_buffers;

	SkmerSaver(std::string outdir) : outdir(outdir)
	{
		// Verify the existance of the out directory
		struct stat sb;

		if (stat(outdir.c_str(), &sb) != 0)
		{
			cerr << "Invalid output folder for " << outdir << endl;
			cerr << "Maybe the folder does not exist ?" << endl;
			exit(1);
		}

		if (outdir[outdir.length() - 1] != '/')
			this->outdir = outdir + "/";
	};

	~SkmerSaver()
	{
		for (auto & pair : this->file_buffers)
		{
			const uint64_t minimizer = pair.first;
			ofstream mini_file(this->outdir + to_string(minimizer) + ".txt");
			string val(pair.second.str());
			mini_file.write(val.c_str(), val.length());
			mini_file.close();
		}
	}

	// Uses the circular buffer from the inputs to save a superkmer in the right buffer.
	void save_skmer (const uint64_t minimizer, const char * buffer, const size_t buff_size, size_t start, size_t stop)
	{
		if (stop == start)
			return;

		// cout << "save " << start << " " << stop << endl;

		// Create a new string stream if not already present
		if (this->file_buffers.find(minimizer) == this->file_buffers.end())
			this->file_buffers[minimizer] = stringstream();

		// Construct the superkmer in the right buffer
		for (size_t i=start ; i<=stop ; i++)
			this->file_buffers[minimizer] << buffer[i % buff_size];
		this->file_buffers[minimizer] << endl;
	};
};
*/

void splitIntoFile(std::string outfile, std::size_t id, const std::size_t k,
                 const std::size_t m, const std::size_t fifo_size, Kmer **fifo,
                 sem_t *empty, sem_t *full) 
{

	std::this_thread::sleep_for(std::chrono::milliseconds(400*(id+1)));

	outfile[outfile.length()-5] = std::to_string(id)[0];

    std::ofstream outf(outfile);

    outf << "bonjour" << "\n";
    
    std::size_t fifo_counter = 0;
    Kmer* sk = new Kmer(2*k-m, false);

    std::cout << "[thread " << id << " start]" << std::endl;

    while(true){
		std::this_thread::sleep_for(std::chrono::milliseconds(400*(id+1)));
        sem_wait(full);

		ssize_t truc = id*fifo_size + fifo_counter;

        sk = fifo[truc];
        fifo_counter = (fifo_counter+1)%fifo_size;

		std::cout << "thread " << id << " in " << (truc) << " recieved " << *sk << std::endl;

        if((*sk).len == 0){ // sent if fasta file is finished
            std::cout<<"!c!"<<std::endl;
            std::cout << "[thread " << id << " over]\n";
            outf.close();
            delete sk;
            return;
        } else {
            outf << sk << "\n";
        }
        sem_post(empty);
    }
}
