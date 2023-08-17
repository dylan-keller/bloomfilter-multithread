#include "../include/KmerAnswer.hpp"

KmerAnswer::KmerAnswer() : bv{}, pos{0}, outbv_nb{0} {}

KmerAnswer::KmerAnswer(std::size_t pos_, std::size_t outbv_nb_) : bv{}, pos{pos_}, outbv_nb{outbv_nb_} {}

void KmerAnswer::fillOutputVector(bm::bvector<>* outbv){
    // outbv->combine_operation_or(bv, pos); // NOTE: this is probably faster, but isn't the right function
    // tried different functions given by BitMagic but no one seems to fit exactly. 
    // for now, let's just do the simple approach :
    for (bm::id_t i=0; i<bv.size(); i++) {
        if (bv.test(i)) outbv->set(pos + i);
    }
}

void KmerAnswer::fillOutputVector(bm::bvector<>* outbv_first, bm::bvector<>* outbv_second, std::size_t bits_first){
    bm::id_t i=0;
    for (i; i<bits_first; i++) { 
        if (bv.test(i)) outbv_first->set(pos + i);
    }
    for (i; i<bv.size(); i++){
        if (bv.test(i)) outbv_second->set(i);
    }
}