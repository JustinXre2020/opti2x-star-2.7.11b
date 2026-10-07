#ifndef H_sjdbInsertSA
#define H_sjdbInsertSA

#include "IncludeDefine.h"
#include "PackedArray.h"
#include <vector>

struct SAinsertReport {
    bool block;     //no entry changed, the entries were moved as a block
    uint64 nRounds; //the work was shared among threads in this many rounds, 0 if one thread did it all
};

// Writes into SA2 the suffix array of the genome with junctions inserted. SA holds the nSA1 entries of the old genome
// (nGenome1 bases, its nSJ1 junctions from sjGstart on, sjdbLength each) and may lie in SA2's buffer past its start.
// Old-junction positions move to their new junctions (oldSJind), minus-strand positions shift by the sjNew new
// junctions, and the nInd new suffixes in indArray (pairs of SA insertion point and position in the 2*nGsj junction
// sequences, sorted, terminated by -999) go before the entries they precede. Up to nThreads threads share the work when
// each of them gets at least minPerThread entries in every round.
SAinsertReport sjdbInsertSAindices(PackedArray &SA, PackedArray &SA2, uint64 nSA1, uint64 nGenome1, uint64 nGenome, uint64 sjGstart,
                                   uint64 sjdbLength, uint64 nSJ1, const uint32 *oldSJind, uint64 sjNew, uint64 nGsj,
                                   const uint64 *indArray, uint64 nInd, char GstrandBit, uint64 nThreads,
                                   uint64 minPerThread=1LLU<<16);

// The rounds in which threads can rebuild the entries: the output entries from rounds[r] to rounds[r+1] are all written
// before the first byte of SA the round reads (SA starts gapBytes past SA2), and rounds[r] is a multiple of 64 for r>0,
// so that the threads of a round split it at whole 64-bit words. Empty if some round would be shorter than minRound.
std::vector<uint64> sjdbInsertSArounds(uint64 gapBytes, uint64 wordLength, uint64 nSA1, const uint64 *indArray, uint64 nInd,
                                       uint64 minRound);

#endif
