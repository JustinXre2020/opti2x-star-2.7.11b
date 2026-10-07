#include "sjdbInsertSA.h"
#include <algorithm>
#include <cstring>

//true if the insertion leaves every old entry where and as it was
static bool entriesUnchanged(uint64 nGenome1, uint64 nGenome, uint64 nSJ1, const uint32 *oldSJind, uint64 sjNew, uint64 nInd) {
    if (nInd>0 || sjNew>0 || nGenome!=nGenome1)
        return false;
    for (uint64 sj1=0; sj1<nSJ1; sj1++) {
        if (oldSJind[sj1]!=sj1)
            return false;
    };
    return true;
};

//writes w-bit entries one after another from the start of a packed array, a whole 64-bit word at a time; the bits past
//the last entry keep their values
struct PackedWriter {
    char *out;
    uint64 w, acc, nBits;
    PackedWriter(char *outIn, uint64 wIn) : out(outIn), w(wIn), acc(0), nBits(0) {};
    inline void put(uint64 x) {//x must fit in w bits
        acc |= x << nBits;
        nBits += w;
        if (nBits>=64) {
            memcpy(out, &acc, 8);
            out += 8;
            nBits -= 64;
            acc = nBits>0 ? x >> (w-nBits) : 0;
        };
    };
    void finish() {//the whole bytes of the last word, then the bits of its last byte
        uint64 nB=nBits/8;
        memcpy(out, &acc, nB);
        if (nBits%8>0) {
            unsigned char keep=(unsigned char) (0xFF << (nBits%8));
            out[nB]=(char) (((unsigned char) out[nB] & keep) | (unsigned char) (acc >> (8*nB)));
        };
    };
};

//the same entries as the upstream loop below, read and written in one pass. Every entry is written behind the entries
//still to be read as long as SA starts far enough past SA2: the loop then reads only old entries, as the upstream one does
static bool readsStayAhead(PackedArray &SA, PackedArray &SA2, uint64 nInd) {
    return SA.charArray>=SA2.charArray && (uint64) (SA.charArray-SA2.charArray) >= nInd*SA.wordLength/8+16;
};

//new entry k goes to output entry min(its insertion point, nSA1)+k; this counts the new entries before output entry J
static uint64 newEntriesBefore(uint64 J, uint64 nSA1, const uint64 *indArray, uint64 nInd) {
    uint64 lo=0, hi=nInd;
    while (lo<hi) {
        uint64 k=(lo+hi)/2;
        if (std::min(indArray[2*k], nSA1)+k < J)
            lo=k+1;
        else
            hi=k;
    };
    return lo;
};

std::vector<uint64> sjdbInsertSArounds(uint64 gapBytes, uint64 wordLength, uint64 nSA1, const uint64 *indArray, uint64 nInd,
                                       uint64 minRound) {
    uint64 nSA2=nSA1+nInd;
    std::vector<uint64> rounds={0};
    for (uint64 J=0; J<nSA2; ) {
        uint64 isa=J-newEntriesBefore(J, nSA1, indArray, nInd);
        uint64 end=(gapBytes+isa*wordLength/8)*8/wordLength;//output entries before this one end before the round's first byte
        if (end>=nSA2)
            end=nSA2;
        else
            end-=end%64;
        if (end<=J || (end<nSA2 && end-J<minRound))
            return {};
        rounds.push_back(end);
        J=end;
    };
    return rounds;
};

//writes output entries J0 to J1; J0 must be a multiple of 64 (or 0), and so must J1 unless it is the last entry
static void insertStreamed(PackedArray &SA, PackedArray &SA2, uint64 nSA1, uint64 nGenome1, uint64 nGenome, uint64 sjGstart,
                           uint64 sjdbLength, const uint32 *oldSJind, uint64 sjNew, uint64 nGsj, const uint64 *indArray,
                           uint64 nInd, char GstrandBit, uint64 J0, uint64 J1) {
    const uint64 w=SA.wordLength, entryMask=(~0LLU) >> (64-w);
    const uint64 N2bit=1LLU << GstrandBit, strandMask=~N2bit, nGsjNew=sjNew*sjdbLength;
    const char *in=SA.charArray;
    PackedWriter out(SA2.charArray+J0*w/8, w);

    auto newEntry = [&](uint64 ind1) -> uint64 {
        return ind1<nGsj ? ind1+sjGstart : (ind1-nGsj) | N2bit;
    };

    uint64 isj=newEntriesBefore(J0, nSA1, indArray, nInd), isa=J0-isj;
    for (uint64 J=J0; J<J1; ) {
        uint64 next= isj<nInd ? std::min(indArray[isj*2], nSA1) : nSA1;
        if (isj<nInd && next==isa) {//insert sj index before the existing index, or past the last one
            out.put(newEntry(indArray[isj*2+1]));
            ++isj; ++J;
            continue;
        };
        uint64 isaEnd=std::min(next, isa+(J1-J));
        J+=isaEnd-isa;
        for (uint64 b=isa*w; isa<isaEnd; isa++, b+=w) {
            uint64 ind1;
            memcpy(&ind1, in+b/8, 8);
            ind1=(ind1 >> (b%8)) & entryMask;
            //the strands of neighboring suffixes are random, so they are told apart with masks rather than branches
            uint64 minusMask=0-(ind1 >> GstrandBit);
            uint64 pos=ind1 ^ (minusMask & (ind1 ^ (nGenome1-(ind1 & strandMask))));//position on the strand of the suffix
            if (__builtin_expect(pos>=sjGstart, 0)) {//this index was an old sj, may need to shift it
                uint64 sj1=(pos-sjGstart)/sjdbLength;//old junction index
                pos+=(oldSJind[sj1]-sj1)*sjdbLength;
                ind1= minusMask ? (nGenome-pos) | N2bit : pos;
            } else {
                ind1+=minusMask & nGsjNew; //reverse complementary indices are all shifted by the length of junctions
            };
            out.put(ind1);
        };
    };
    out.finish();
};

//rounds of output entries, each shared among the threads at multiples of 64 entries; returns the number of rounds, or 0
//if the room between SA2 and SA is too short for rounds of minPerThread entries per thread
static uint64 insertShared(PackedArray &SA, PackedArray &SA2, uint64 nSA1, uint64 nGenome1, uint64 nGenome, uint64 sjGstart,
                           uint64 sjdbLength, const uint32 *oldSJind, uint64 sjNew, uint64 nGsj, const uint64 *indArray,
                           uint64 nInd, char GstrandBit, uint64 nThreads, uint64 minPerThread) {
    if (nThreads<2 || SA.charArray<SA2.charArray)
        return 0;
    std::vector<uint64> rounds=sjdbInsertSArounds(SA.charArray-SA2.charArray, SA.wordLength, nSA1, indArray, nInd, minPerThread*nThreads);
    if (rounds.empty())
        return 0;
    for (uint64 r=0; r+1<rounds.size(); r++) {
        uint64 J0=rounds[r], span=rounds[r+1]-J0;
        #pragma omp parallel for num_threads(nThreads) schedule(static,1)
        for (uint64 t=0; t<nThreads; t++) {
            uint64 a=J0+span*t/nThreads/64*64, b= t+1==nThreads ? J0+span : J0+span*(t+1)/nThreads/64*64;
            if (a<b)
                insertStreamed(SA, SA2, nSA1, nGenome1, nGenome, sjGstart, sjdbLength, oldSJind, sjNew, nGsj, indArray, nInd, GstrandBit, a, b);
        };
    };
    return rounds.size()-1;
};

//moves nBytes from SA to SA2 in rounds of the room between them: a round writes only where earlier rounds read, so its
//threads copy their parts at once; returns the number of rounds, or 0 if the bytes were moved by one thread
static uint64 moveShared(PackedArray &SA, PackedArray &SA2, uint64 nBytes, uint64 nThreads, uint64 minPerThread) {
    uint64 gap= SA.charArray>=SA2.charArray ? SA.charArray-SA2.charArray : 0;
    if (nThreads<2 || gap*8/SA.wordLength < minPerThread*nThreads) {
        memmove(SA2.charArray, SA.charArray, nBytes);
        return 0;
    };
    uint64 nRounds=(nBytes+gap-1)/gap;
    for (uint64 r=0; r<nRounds; r++) {
        uint64 b0=r*gap, span=std::min(nBytes-b0, gap);
        #pragma omp parallel for num_threads(nThreads) schedule(static,1)
        for (uint64 t=0; t<nThreads; t++) {
            uint64 a=b0+span*t/nThreads, b=b0+span*(t+1)/nThreads;
            memcpy(SA2.charArray+a, SA.charArray+a, b-a);
        };
    };
    return nRounds;
};

SAinsertReport sjdbInsertSAindices(PackedArray &SA, PackedArray &SA2, uint64 nSA1, uint64 nGenome1, uint64 nGenome, uint64 sjGstart,
                                   uint64 sjdbLength, uint64 nSJ1, const uint32 *oldSJind, uint64 sjNew, uint64 nGsj,
                                   const uint64 *indArray, uint64 nInd, char GstrandBit, uint64 nThreads, uint64 minPerThread) {

    if (entriesUnchanged(nGenome1, nGenome, nSJ1, oldSJind, sjNew, nInd)) {
        //move the bytes up to the one the last entry starts in, all of whose bits before the last entry belong to entries;
        //the last entry is written on its own, so the bits past it keep their values, as with entry-by-entry writes
        uint64 last=SA[nSA1-1];
        uint64 nRounds=moveShared(SA, SA2, (nSA1-1)*SA.wordLength/8+1, nThreads, minPerThread);
        SA2.writePacked(nSA1-1,last);
        return {true, nRounds};
    };

    uint64 nRounds=insertShared(SA, SA2, nSA1, nGenome1, nGenome, sjGstart, sjdbLength, oldSJind, sjNew, nGsj, indArray, nInd,
                                GstrandBit, nThreads, minPerThread);
    if (nRounds>0)
        return {false, nRounds};

    if (readsStayAhead(SA, SA2, nInd)) {
        insertStreamed(SA, SA2, nSA1, nGenome1, nGenome, sjGstart, sjdbLength, oldSJind, sjNew, nGsj, indArray, nInd, GstrandBit,
                       0, nSA1+nInd);
        return {false, 0};
    };

    uint nGsjNew=sjNew*sjdbLength; //this is the actual number of bytes added to the genome, while nGsj is the total size of all junctions

    uint N2bit= 1LLU << GstrandBit;
    uint strandMask=~N2bit;

    uint isj=0, isa2=0;
    for (uint isa=0;isa<nSA1;isa++) {
        while (isa==indArray[isj*2]) {//insert sj index before the existing index
            uint ind1=indArray[isj*2+1];
            if (ind1<nGsj) {
                ind1+=sjGstart;
            } else {//reverse strand
                ind1=(ind1-nGsj) | N2bit;
            };
            SA2.writePacked(isa2,ind1);
            ++isa2; ++isj;
        };

        uint ind1=SA[isa];
        if ( (ind1 & N2bit)>0 )
        {//- strand
            uint ind1s = nGenome1 - (ind1 & strandMask);
            if (ind1s>=sjGstart)
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1s-sjGstart)/sjdbLength;//old junction index
                ind1s += (oldSJind[sj1]-sj1)*sjdbLength;
                ind1 = (nGenome - ind1s) | N2bit;
            } else
            {
                ind1+=nGsjNew; //reverse complementary indices are all shifted by the length of junctions
            };
        } else
        {//+ strand
            if (ind1>=sjGstart)
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1-sjGstart)/sjdbLength;//old junction index
                ind1 += (oldSJind[sj1]-sj1)*sjdbLength;
            };
        };

        SA2.writePacked(isa2,ind1);
        ++isa2;
    };

    for (;isj<nInd;isj++) {//insert last new indices after the last old index
        uint ind1=indArray[isj*2+1];
        if (ind1<nGsj) {
            ind1+=sjGstart;
        } else {//reverse strand
            ind1=(ind1-nGsj) | N2bit;
        };
        SA2.writePacked(isa2,ind1);
        ++isa2;
    };
    return {false, 0};
};
