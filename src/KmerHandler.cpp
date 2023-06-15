#include "../include/KmerHandler.hpp"

KmerHandler::KmerHandler(std::size_t k_) :
    k(k_), current_kmer(0), mask((1ull << (2 * (k-1))) - 1)
{};

KmerHandler::~KmerHandler()
{}

__uint128_t KmerHandler::next_Kmer(char nucleotide){
    // Remove the leftmost nucleotide
        // e.g. 11111111 becomes 00111111
    this->current_kmer &= this->mask;

    // Shift the kmer to the left creating a "hole" on the right
        // e.g. 00111111 becomes 11111100
    this->current_kmer <<= 2;

    // Add the new nucleotide
        // e.g. 11111100 becomes 111111XX, XX being the bits of nucleotid
    switch (nucleotide) {
    //case 'A': // not sure if removing this speeds up or if the compiler knows
    //case 'a':
    //    this->current_kmer += 0;
    //    break;
    case 'C':
    case 'c':
        this->current_kmer += 1;
        break;
    case 'G':
    case 'g':
        this->current_kmer += 2;
        break;
    case 'T':
    case 't':
        this->current_kmer += 3;
        break;
    }

    return this->current_kmer;
}

std::string uint_To_Kmer (__uint128_t val, std::size_t k)
{
	std::stringstream ss;
	for (std::size_t i=0 ; i<k ; i++)
	{
		char c;
		switch (val & 0b11) {
		case 0: c = 'A'; break;
		case 1: c = 'C'; break;
		case 2: c = 'G'; break;
		case 3: c = 'T'; break;
		}
		ss << c;
		val >>= 2;
	}

	std::string kmer = ss.str();
	reverse(kmer.begin(), kmer.end());
	return kmer;
}