#ifndef FASTAREADER_HPP
#define FASTAREADER_HPP

#include <string>
#include <fstream>
#include <iostream>

class FastaReader{
  public:
    std::ifstream fs;
    
    FastaReader(std::string filename);
    ~FastaReader();

    bool has_next();

    char next_char();
};

#endif