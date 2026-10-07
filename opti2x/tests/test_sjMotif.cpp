// SjMotifTable must classify intron-end bases exactly like the if-chain in upstream stitchAlignToTranscript.cpp.
#include "sjMotif.h"
#include <cstdio>
#include <random>

// Verbatim motif chain from upstream stitchAlignToTranscript.cpp (deletion branch); d1 d2 = G[gAend+jR1+1..2],
// a1 a2 = G[gBstart1+jR1-1..0].
static int upstreamMotif(char d1, char d2, char a1, char a2) {
    if ( d1==2 && d2==3 && a1==0 && a2==2 ) {//GTAG
        return 1;
    } else if ( d1==1 && d2==3 && a1==0 && a2==1 ) {//CTAC
        return 2;
    } else if ( d1==2 && d2==1 && a1==0 && a2==2 ) {//GCAG
        return 3;
    } else if ( d1==1 && d2==3 && a1==2 && a2==1 ) {//CTGC
        return 4;
    } else if ( d1==0 && d2==3 && a1==0 && a2==1 ) {//ATAC
        return 5;
    } else if ( d1==2 && d2==3 && a1==0 && a2==3 ) {//GTAT
        return 6;
    } else {
        return 0;
    };
}

static int candidateMotif(char d1, char d2, char a1, char a2) {
    return sjMotifTable(SjMotifTable::base(d1), SjMotifTable::base(d2), SjMotifTable::base(a1), SjMotifTable::base(a2));
}

static int nFail=0;
static void check(char d1, char d2, char a1, char a2) {
    int u=upstreamMotif(d1,d2,a1,a2), c=candidateMotif(d1,d2,a1,a2);
    if (u!=c && nFail++<5)
        printf("  motif(%d,%d,%d,%d): upstream %d, candidate %d\n", d1, d2, a1, a2, u, c);
}

static bool testEveryGenomeCodeCombination() {
    //genome codes: ACGT 0-3, N 4, spacer 5, fragment spacer 11, and arbitrary bytes incl. negative chars
    const char codes[]={0,1,2,3,4,5,11,127,-1,-128};
    nFail=0;
    for (char d1: codes) for (char d2: codes) for (char a1: codes) for (char a2: codes)
        check(d1,d2,a1,a2);
    return nFail==0;
}

static bool testRandomBytes() {
    std::mt19937 rng(20261007);
    nFail=0;
    for (int i=0; i<2000000; i++) {
        uint32_t x=rng();
        check((char)x, (char)(x>>8), (char)(x>>16), (char)(x>>24));
    }
    return nFail==0;
}

int main() {
    int rc=0;
    if (!testEveryGenomeCodeCombination()) {printf("FAIL testEveryGenomeCodeCombination (%d)\n", nFail); rc=1;}
    if (!testRandomBytes()) {printf("FAIL testRandomBytes (%d)\n", nFail); rc=1;}
    return rc;
}
