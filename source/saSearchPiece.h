#ifndef H_saSearchPiece
#define H_saSearchPiece

#include "IncludeDefine.h"
#include "Genome.h"

// Upper end of the suffix array range of the present SAindex prefix ind1 of length Lind; good=false when it is only a
// safe bound, not the range end
uint saiUpperBound(Genome &mapGen, uint Lind, uint ind1, bool& good);

// Searches the suffix array for the longest genome match of one read piece: Read1[0] from pieceStart forward if dirR,
// else backward on the complement. Returns the number of matches, sets maxL and the SA range in indStartEnd.
uint searchPiece(Genome &mapGen, char** Read1, uint pieceStart, uint pieceLength, bool dirR, uint& maxL, uint* indStartEnd);

#endif
