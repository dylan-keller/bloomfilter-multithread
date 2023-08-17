#include "../include/SkmerExtractor.hpp"

// TODO : remove "truc" name

void extractSkmers(std::string filename, const std::size_t k, const std::size_t m,
                const std::size_t q, const std::size_t fifo_size, Kmer** fifos,
                sem_t* emptys, sem_t* fulls, std::atomic<std::size_t>* end_increment){
    
    // ------------------------------ Variables ------------------------------ 

    FastaReader fr(filename);

    // Our k-sized window as we read the text, a.k.a. our "current" k-mer
    std::string kmer_cur;
    kmer_cur.resize(k);

    // The forward (resp. reverse-strand) hash values of
    // every m-mer in the k-mer, for quick minimizer find
    //     reminder : a k-mer contains (k-m+1) m-mers = possible minimizers
    uint64_t fhvalues[k-m+1];
    uint64_t rhvalues[k-m+1];

    // The forward (resp. reverse-strand) hash values
    uint64_t hVal=0, fhVal=0, rhVal=0; 

    // Tells us our kmer is reverse complement or not
    bool isRevComp;

    // The minimizer hash values :
    uint64_t fhpos; // position (in the array) of the forward minimizer
    uint64_t rhpos; // position (in the array) of the reverse-strand minimizer
    uint64_t hpos;  // position (in the array) of the canonical minimizer
    uint64_t fhmin; // hash value of the forward minimizer
    uint64_t rhmin; // hash value of the reverse-strand minimizer 
    uint64_t hmin;  // hash value of the canonical minimizer 

    // Loop counter (used for our circular arrays fhvalues and rhvalues) and normal counter
    uint16_t counter;
    uint64_t actual_counter = 0;
    
    // Threads, each associated to a super-k-mer fifo
    std::vector<std::thread> thread_fifos;

    // Flag if the super-k-mer has ended and we need to create a new one
    bool new_skmer_flag = false;

    // Last read char
    char c;

    // Used to know the fifo in which to add the super-k-mer
    std::size_t fifo_nb;

    // Used to know where we are in our circular arrays
    std::size_t fifo_counter[q];
    for (std::size_t i=0; i<q; i++) fifo_counter[i]=0;

    // ------------------------------ Program ------------------------------

    do {
        // To start, we need to read all the k first characters 
        for (std::size_t i=0; i<k; i++){
            kmer_cur[i] = fr.next_char();
        }

        hVal = NTC64(kmer_cur.c_str(), m, fhVal, rhVal);
        // std::cout << "[" << fhVal << "_" << rhVal << "]\n";
        fhvalues[0] = fhVal;
        rhvalues[0] = rhVal;

        // Let's compute the hash values of all m-mers in our first k-mer using ntHash
        for (std::size_t i=0; i<k-m; i++){
            NTC64(kmer_cur[i], kmer_cur[i+m], m, fhVal, rhVal);
            fhvalues[i+1] = fhVal;
            rhvalues[i+1] = rhVal;
            // std::cout << "{" << i << "_" << fhVal << "_" << rhVal << "}\n";
        }

        // Among these m-mers, let's find the minimizer (our first minimizer)
        fhpos = std::min_element(fhvalues, fhvalues+k-m+1) - fhvalues;
        rhpos = std::min_element(rhvalues, rhvalues+k-m+1) - rhvalues;
        fhmin = fhvalues[fhpos];
        rhmin = rhvalues[rhpos];
        
        if (fhmin < rhmin){ 
            hmin = fhmin;
            hpos = fhpos;
            isRevComp = false;
        } else { // note : if they are equal, we arbitrarily pick the forward one
            hmin = rhmin;
            hpos = rhpos;
            isRevComp = true;
        } 

        // We create the first super-k-mer with our current k-mer
        Kmer* sk = new Kmer(2*k-m, actual_counter, isRevComp, kmer_cur);

        // And we start reading the rest of the file
        c = fr.next_char();
        counter = 0;
        actual_counter++;

        // for(int ii=0; ii<1100; ii++){ // TEST
        while(c != '\0'){ // \0 should be returned at the end of a sequence

            // Get the hash value of the new m-mer (rightmost)
            NTC64(kmer_cur[k-m], c, m, fhVal, rhVal);

            // Get the next k-mer (rotate the std::string once leftwise, and replace last character)
            rotate(kmer_cur.begin(), kmer_cur.begin()+1, kmer_cur.end());
            kmer_cur[k-1] = c;

            // if (actual_counter<5){
            //     std::cout << '{' <<  kmer_cur << '_' << fhVal << '_' << rhVal << "}\n";
            // }

            // Place the new hash values in their respective arrays
            // Thanks to 'counter', we can use the array as a circular array
            fhvalues[counter] = fhVal;
            rhvalues[counter] = rhVal;

            /*
             * Now that we read a new character and therefore advanced the k-sized window,
             * we have many things to check, to see if it's the end of the current super-k-mer or not.
             * The possibilites are :
             * 
             * (1) the previous minimizer fell out of the k-sized window (if was the leftmost minimizer)
             * (2) the new m-mer (rightmost) is a better minimizer than the previous one
             * (3) the new m-mer (rightmost) is not better and we keep the previous one
             * ((4)) we reached the end of the file
             */

            // Let's check if the prev minimizer is still part of the new k-mer.
            if (counter != hpos){ // If it is,
                // We need to check if the new m-mer (rightmost) is a better minimizer.
                if (fhvalues[counter] < hmin) {
                    hmin = fhvalues[counter];
                    hpos = counter;
                    isRevComp = false;
                    new_skmer_flag = true;
                } else if (rhvalues[counter] < hmin) {
                    hmin = rhvalues[counter];
                    hpos = counter;
                    isRevComp = true;
                    new_skmer_flag = true;
                } else {
                    sk->addNucl(c);
                }
            } 
            else { // If the minimizer fell out of the window, we end the super-k-mer.
                new_skmer_flag = true;

                fhpos = std::min_element(fhvalues, fhvalues+k-m+1) - fhvalues;
                rhpos = std::min_element(rhvalues, rhvalues+k-m+1) - rhvalues;
                fhmin = fhvalues[fhpos];
                rhmin = rhvalues[rhpos];

                if (fhmin < rhmin){ 
                    hmin = fhmin;
                    hpos = fhpos;
                    isRevComp = false;
                } else { // note : if they are equal, we arbitrarily pick the forward one
                    hmin = rhmin;
                    hpos = rhpos;
                    isRevComp = true;
                } 
            }

            if (new_skmer_flag){
                // We need to add the super-k-mer to its correct fifo.
                fifo_nb = hmin%q;
                // Wait for an empty spot i, the correct fifo
                sem_wait(&emptys[fifo_nb]);
                // add the super-k-mer in the correct spot
                ssize_t truc = fifo_nb*fifo_size + fifo_counter[fifo_nb];
                fifos[truc] = sk;
                // update the counter
                fifo_counter[fifo_nb] = (fifo_counter[fifo_nb]+1)%fifo_size;
                // Notifies that a spot has just been filled
                sem_post(&fulls[fifo_nb]);
                // We can now create the new super-k-mer, starting at the current k-mer.
                sk = new Kmer(2*k-m, actual_counter, isRevComp, kmer_cur);
            }

            c = fr.next_char();
            counter = (counter+1)%(k-m+1);
            actual_counter++;
            new_skmer_flag = false;
        }

        // when the sequence (or file) is over, we need to send the final super-k-mer
        
        // cf. lines just above if clarification needed ("if(new_skmer_flag)...")
        fifo_nb = hmin%q;
        sem_wait(&emptys[fifo_nb]);
        ssize_t truc = fifo_nb*fifo_size + fifo_counter[fifo_nb];
        //std::cout << "put to " << fifo_nb << " in " << truc << ": " << *sk << std::endl;
        fifos[truc] = sk;
        fifo_counter[fifo_nb] = (fifo_counter[fifo_nb]+1)%fifo_size;
        sem_post(&fulls[fifo_nb]);
    
    } while (false); // TEST ; for now we don't want to loop over the whole file
    //} while (fr.has_next());

    std::cout << "(((extraction done)))" << std::endl;

    for (std::size_t i=0; i<q; i++){
        Kmer* kmer_ender = new Kmer(1, false);
        //std::cout << "sending kill signal to " << i << std::endl;
        sem_wait(&emptys[i]);
        fifos[i*fifo_size + fifo_counter[i]] = kmer_ender;
        sem_post(&fulls[i]);
    }

    if (end_increment != nullptr){
        end_increment[0]++;
    }

    std::cout << "[extractor thread over]" << std::endl;
}