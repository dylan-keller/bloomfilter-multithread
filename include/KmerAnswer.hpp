#ifndef KMERANSWER_HPP
#define KMERANSWER_HPP

#include <iostream>

#include "/mnt/c/Users/dylan/Documents/M1/stage/projet2/bloomfilter-multithread/external/bitmagic/src/bm.h"

class KmerAnswer{
  public :
    bm::bvector<> bv;
    bm::id_t pos;
    std::size_t outbv_nb;

    KmerAnswer();
    KmerAnswer(std::size_t pos_, std::size_t outbv_nb_);

    // Fills the given bitvector at position 'pos' with the bits of 'bv'.
    void fillOutputVector(bm::bvector<>* outbv);

    // Used if the bits to fill are overlapping two bitvectors.
    //   (example: bv.size = 10, pos = 995, and outbv.size = 1000. not enough room for all bits)
    // Fills the 'outbv_first' at position 'pos' with the 'bits_first' first bits of 'bv', 
    // then fills the remaining bits at position 0 of 'outbv_second'
    void fillOutputVector(bm::bvector<>* outbv_first, bm::bvector<>* outbv_second, std::size_t bits_first);
};

#endif