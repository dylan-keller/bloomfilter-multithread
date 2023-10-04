#include "../include/KmerAnswer.hpp"

KmerAnswer::KmerAnswer() : bv{}, pos{0}, outbv_nb{0} {}

KmerAnswer::KmerAnswer(std::size_t pos_, std::size_t outbv_nb_) : bv{}, pos{pos_}, outbv_nb{outbv_nb_} {}

void KmerAnswer::fillOutputVector(bm::bvector<>* outbv){
    // outbv->combine_operation_or(bv, pos); // NOTE: this is probably faster, but isn't the right function
    // tried different functions given by BitMagic but no one seems to fit exactly. 
    // for now, let's just do the simple approach :
    for (bm::id_t i=0; i<bv.size(); i++) if (bv.test(i)) outbv->set(pos + i);
}

void KmerAnswer::fillOutputVector(bm::bvector<>* outbv_first, bm::bvector<>* outbv_second, std::size_t outbv_size){
    bm::id_t i=0;

    if (pos+bv.size() >= outbv_size){ // if this KmerAnswer overlaps two bitvectors :
        // this code is ugly, but this operation is very rare
        for (bm::id_t i=0; i<bv.size(); i++){
            if (bv.test(i)){
                if (pos+i < outbv_size) outbv_first->set(pos + i);
                else outbv_second->set(pos+i-outbv_size);
            }
        }
    } else { // otherwise, simple fill
        for (bm::id_t i=0; i<bv.size(); i++) if (bv.test(i)) outbv_first->set(pos + i);
    }
}

void KmerAnswer::fillOutputVector(bm::bvector<>* outbv_first, bm::bvector<>* outbv_second, std::size_t outbv_size,
                                  std::size_t* counter_first, std::size_t* counter_second){
    bm::id_t i=0;

    if (pos+bv.size() >= outbv_size){ // if this KmerAnswer overlaps two bitvectors :
        std::size_t first_bits = outbv_size - pos;
        for(; i<first_bits; i++){
            if (bv.test(i)) outbv_first->set(pos + i);
        }
        for(; i<bv.size(); i++){
            if (bv.test(i)) outbv_second->set(i-first_bits); // write the rest in second bitvector starting at 0
        }
        *counter_first += first_bits;
        *counter_second += bv.size() - first_bits;
    } else { // otherwise, simple fill
        for (; i<bv.size(); i++) if (bv.test(i)) outbv_first->set(pos + i);
        *counter_first += bv.size();
    }
}