#include "ReadAlign.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"
#include "saSearchPiece.h"

uint ReadAlign::maxMappableLength2strands(uint pieceStartIn, uint pieceLengthIn, uint iDir, uint /*iSA1*/, uint /*iSA2*/, uint& maxLbest, uint iFrag) {
    //returns number of mappings, maxMappedLength=mapped length
    uint Nrep=0, indStartEnd[2], maxL;

    uint NrepAll[P.pGe.gSAsparseD], indStartEndAll[P.pGe.gSAsparseD][2], maxLall[P.pGe.gSAsparseD];
    maxLbest=0;

    bool dirR = iDir==0;

    // defaults:  (from genomeParameters.txt)
    // gSAsparseD = 1
    // gSAindexNbases = 14

    for (uint iDist=0; iDist<min(pieceLengthIn,P.pGe.gSAsparseD); iDist++) {//cycle through different distances
        uint pieceLength=pieceLengthIn-iDist;
        uint pieceStart = dirR ? pieceStartIn+iDist : pieceStartIn-iDist;
        Nrep = searchPiece(mapGen, Read1, pieceStart, pieceLength, dirR, maxL, indStartEnd);

        if (maxL+iDist > maxLbest) {//this idist is better
            maxLbest=maxL+iDist;
        };
        NrepAll[iDist]=Nrep;
        indStartEndAll[iDist][0]=indStartEnd[0];
        indStartEndAll[iDist][1]=indStartEnd[1];
        maxLall[iDist]=maxL;
    };

    for (uint iDist=0; iDist<min(pieceLengthIn,P.pGe.gSAsparseD); iDist++) {//cycle through different distances, store the ones with largest maxL
        if ( (maxLall[iDist]+iDist) == maxLbest) {
            storeAligns(iDir, (dirR ? pieceStartIn+iDist : pieceStartIn-iDist), NrepAll[iDist], maxLall[iDist], indStartEndAll[iDist], iFrag);
        };
    };
    return Nrep;
};
