#ifndef H_WinBinReset
#define H_WinBinReset

#include "IncludeDefine.h"

// Restores the winBin table to empty (uintWinBinMax) between reads by resetting only the bin intervals
// written since the previous reset; falls back to resetting the whole table if there were too many.
// For a human-size genome the whole table is ~190KB, which dominated stitchPieces for short reads.
class WinBinReset {
public:
    static const uint nMax=512;

    void init(uintWinBin **winBinIn, uint winBinNin) {
        winBin=winBinIn;
        winBinN=winBinNin;
        n=0;
    };

    void mark(uint str, uint b1, uint b2) {//bins b1..b2 (b1<=b2) of strand str were written
        if (n>nMax)
            return;
        if (n>0) {//merge with the previous interval if it touches this one
            uint *d=dirty[n-1];
            if (d[0]==str && b1<=d[2]+1 && b2+1>=d[1]) {
                d[1]=min(d[1],b1);
                d[2]=max(d[2],b2);
                return;
            };
        };
        if (n==nMax) {
            n++;
            return;
        };
        dirty[n][0]=str;
        dirty[n][1]=b1;
        dirty[n][2]=b2;
        n++;
    };

    void reset() {
        if (n>nMax) {
            memset(winBin[0],255,sizeof(winBin[0][0])*winBinN);
            memset(winBin[1],255,sizeof(winBin[0][0])*winBinN);
        } else {
            for (uint ii=0; ii<n; ii++)
                memset(winBin[dirty[ii][0]]+dirty[ii][1],255,sizeof(winBin[0][0])*(dirty[ii][2]-dirty[ii][1]+1));
        };
        n=0;
    };

private:
    uintWinBin **winBin;
    uint winBinN;
    uint n; //number of recorded intervals; n>nMax means the list overflowed
    uint dirty[nMax][3]; //strand, first bin, last bin
};

#endif
