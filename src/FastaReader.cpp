#include "../include/FastaReader.hpp"

FastaReader::FastaReader(std::string filename){
    this->fs.open(filename);
    if (! this->fs.is_open()) {
        std::cerr << "Impossible to open " << filename << std::endl;
        exit(1);
    }

    // Read the first header
    char c = '\0';
    while (c != '\n') {
        this->fs.get(c);
    }
}

FastaReader::~FastaReader(){
    this->fs.close();
}

bool FastaReader::has_next() {
    return ! this->fs.eof();
}

char FastaReader::next_char() {
    while (true) {
        if (! this->has_next())
            return '\0';

        char c = '\0'; 
        this->fs.get(c);

        switch (c) {
        case 'A':
        case 'a':
        case 'C':
        case 'c':
        case 'G':
        case 'g':
        case 'T':
        case 't':
            return c;
        case '>':
            while (c != '\n' and this->has_next())
                this->fs.get(c);
            return '\0';
        }
    }
}