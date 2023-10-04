# Bloomfilter Multi-thread

This is a C++ implementation of a multi-threaded Bloom filter. It uses two external libraries: ntHash-AVX512 for its rolling hash function, and BitMagic for its efficient bitvectors.

# Usage

Input files :
- one file containing the DNA data that needs to be queried, in FASTA format.
- one file containing the sequences that we want to query, in FASTA format.

Output :

A bitvector denoting the presence (1) or absence (0) of each single k-mer of the input query, in order.

To compile a first time :
```
make
```

To run :
```
// Redaction of this README is in progress.
```

To clean up the repository and delete all compiled files :
```
make clean
```
