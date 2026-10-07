#include "saSearchPiece.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"

uint saiUpperBound(Genome &mapGen, uint Lind, uint ind1, bool& good) {
    good = true;
    uint iSA2;
    if (mapGen.genomeSAindexStart[Lind-1]+ind1+1 < mapGen.genomeSAindexStart[Lind]) {//we are not at the end of the SA
        //an absent prefix holds the SA start of the next present one (nSA past the last), so the bound is as tight as
        //upstream's present-neighbour bound, and the SA search no longer runs to the end of the suffix array
        iSA2 = (mapGen.SAi[mapGen.genomeSAindexStart[Lind-1]+ind1+1] & mapGen.SAiMarkNmask & mapGen.SAiMarkAbsentMask) - 1;
    } else {
        iSA2=mapGen.nSA-1;
        good = false;
    };
    return iSA2;
};

uint searchPiece(Genome &mapGen, char** Read1, uint pieceStart, uint pieceLength, bool dirR, uint& maxL, uint* indStartEnd) {
    Parameters &P = mapGen.P;
    uint Nrep=0, iSA1, iSA2;

    //calculate full index
    uint Lmax=min(mapGen.pGe.gSAindexNbases,pieceLength);
    uint ind1=0;
    if (dirR) {//forward search
        for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
            ind1 <<=2LLU;
            ind1 += ((uint) Read1[0][pieceStart+ii]);
        };
    } else {//reverse search
        for (uint ii=0;ii<Lmax;ii++) {//calculate index TODO: make the index calculation once for the whole read and store it
            ind1 <<=2LLU;
            ind1 += ( 3-((uint) Read1[0][pieceStart-ii]) );
        };
    };

    //find SA boundaries
    uint Lind=Lmax;
    while (Lind>0) {//check the presence of the prefix for Lind
        iSA1=mapGen.SAi[mapGen.genomeSAindexStart[Lind-1]+ind1]; // starting point for suffix array search.
        if ((iSA1 & mapGen.SAiMarkAbsentMaskC) == 0) {//prefix exists
            break;
        } else {//this prefix does not exist, reduce Lind
            --Lind;
            ind1 = ind1 >> 2;
        };
    };

    // define upper bound for suffix array range search.
    bool iSA2good;
    iSA2 = saiUpperBound(mapGen, Lind, ind1, iSA2good);

//#define SA_SEARCH_FULL

#ifdef SA_SEARCH_FULL
    //full search of the array even if the index search gave maxL
    maxL=0;
    Nrep = maxMappableLength(mapGen, Read1, pieceStart, pieceLength, iSA1 & mapGen.SAiMarkNmask, iSA2, dirR, maxL, indStartEnd);
#else
    bool iSA1noN = (iSA1 & mapGen.SAiMarkNmaskC)==0;
    if (Lind < mapGen.pGe.gSAindexNbases && iSA1noN && iSA2good) {//no need for SA search
        // very short seq, already found hits in suffix array w/o having to search the genome for extensions.
        indStartEnd[0]=iSA1;
        indStartEnd[1]=iSA2;
        Nrep=indStartEnd[1]-indStartEnd[0]+1;
        maxL=Lind;
    } else if (iSA1==iSA2 && iSA1noN && iSA2good) {//unique align already, just find maxL
        if ((iSA1 & mapGen.SAiMarkNmaskC)!=0) {
            ostringstream errOut;
            errOut  << "BUG: in ReadAlign::maxMappableLength2strands";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_BUG, P);
        };
        indStartEnd[0]=indStartEnd[1]=iSA1;
        Nrep=1;
        bool comparRes;
        maxL=compareSeqToGenome(mapGen, Read1, pieceStart, pieceLength, Lind, iSA1, dirR, comparRes);
    } else {//need SA search, pieceLength>maxL
        if (iSA2good && iSA1noN) {
            maxL = Lind; //Lind bases were already matched
        } else {
            maxL=0;
        };
        Nrep = maxMappableLength(mapGen, Read1, pieceStart, pieceLength, iSA1 & mapGen.SAiMarkNmask, iSA2, dirR, maxL, indStartEnd);
    };
#endif
    return Nrep;
};
