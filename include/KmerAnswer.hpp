#ifndef KMERANSWER_HPP
#define KMERANSWER_HPP

#include <iostream>

#include "../external/bitmagic/src/bm.h"

class KmerAnswer{
  public :
    bm::bvector<> bv;
    bm::id_t pos;
    std::size_t outbv_nb;

    KmerAnswer();
    KmerAnswer(std::size_t pos_, std::size_t outbv_nb_);

    // Fills the given bitvector at position 'pos' with the bits of 'bv'.
    void fillOutputVector(bm::bvector<>* outbv);

    // Fills the 'outbv_first' bitvector at position 'pos' with the bits from 'bv'.
    // if 'bv' is written near the end of 'outbv_first' and is too big to fit in it,
    // the first bits of 'bv' are written in 'outbv_first' and the rest in 'outbv_second' starting at position 0
    void fillOutputVector(bm::bvector<>* outbv_first, bm::bvector<>* outbv_second, std::size_t outbv_size);

    void fillOutputVector(bm::bvector<>* outbv_first, bm::bvector<>* outbv_second, std::size_t outbv_size,
                          std::size_t* counter_first, std::size_t* counter_second);
};

#endif