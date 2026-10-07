// Differential test: rebuilding the suffix array after junction insertion must leave every byte of the buffer exactly as
// upstream STAR 2.7.11b's sjdbBuildIndex loop (copied below verbatim, minus its commented-out testing code) leaves it,
// with SA lying at the tail of SA2's buffer as genomeLoad lays them out. When no entry changes - no new suffixes, no new
// junctions, the same old junctions in their places - the entries must be moved as a block. Threads share the work in
// rounds whenever the room between SA2 and SA allows; each round's output must lie before the bytes the round reads.
#include "sjdbInsertSA.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

static std::mt19937_64 rng(20261007);
static int nFail=0;
static uint64 nMultiRound=0, nBlockShared=0;//cases in which threads rebuilt the entries in 2+ rounds, or shared a block move

struct GenomeLite {//the Genome fields the upstream loop reads
    uint nSA, nGenome, sjdbLength, nChrReal;
    uint chrStart[2];
    char GstrandBit;
};

namespace ref {
static void insertSA(PackedArray &SA, PackedArray &SA2, GenomeLite &mapGen, GenomeLite &mapGen1, uint32 *oldSJind, uint64 *indArray,
                     uint nInd, uint nGsj, uint64 sjNew) {
    uint nGsjNew=sjNew*mapGen.sjdbLength; //this is the actual number of bytes added to the genome, while nGsj is the total size of all junctions

    uint N2bit= 1LLU << mapGen.GstrandBit;
    uint strandMask=~N2bit;

    uint isj=0, isa2=0;
    for (uint isa=0;isa<mapGen1.nSA;isa++) {
        while (isa==indArray[isj*2]) {//insert sj index before the existing index
            uint ind1=indArray[isj*2+1];
            if (ind1<nGsj) {
                ind1+=mapGen.chrStart[mapGen.nChrReal];
            } else {//reverse strand
                ind1=(ind1-nGsj) | N2bit;
            };
            SA2.writePacked(isa2,ind1);
            ++isa2; ++isj;
        };

        uint ind1=SA[isa];
        if ( (ind1 & N2bit)>0 )
        {//- strand
            uint ind1s = mapGen1.nGenome - (ind1 & strandMask);
            if (ind1s>=mapGen.chrStart[mapGen.nChrReal])
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1s-mapGen.chrStart[mapGen.nChrReal])/mapGen.sjdbLength;//old junction index
                ind1s += (oldSJind[sj1]-sj1)*mapGen.sjdbLength;
                ind1 = (mapGen.nGenome - ind1s) | N2bit;
            } else
            {
                ind1+=nGsjNew; //reverse complementary indices are all shifted by the length of junctions
            };
        } else
        {//+ strand
            if (ind1>=mapGen.chrStart[mapGen.nChrReal])
            {//this index was an old sj, may need to shift it
                uint sj1 = (ind1-mapGen.chrStart[mapGen.nChrReal])/mapGen.sjdbLength;//old junction index
                ind1 += (oldSJind[sj1]-sj1)*mapGen.sjdbLength;
            };
        };

        SA2.writePacked(isa2,ind1);
        ++isa2;
    };

    for (;isj<nInd;isj++) {//insert last new indices after the last old index
        uint ind1=indArray[isj*2+1];
        if (ind1<nGsj) {
            ind1+=mapGen.chrStart[mapGen.nChrReal];
        } else {//reverse strand
            ind1=(ind1-nGsj) | N2bit;
        };
        SA2.writePacked(isa2,ind1);
        ++isa2;
    };
}
}

enum Kind {UNCHANGED, GENERAL, NEW_JUNCTIONS_ONLY, REORDERED, DROPPED};

struct Case {
    Kind kind;
    char GstrandBit;
    uint64 nSA1, sjGstart, sjdbLength, nSJ1, nSJ, sjNew, nInd, extra;
    std::vector<uint32> oldSJind;
    std::vector<uint64> indArray, sa;
};

//old junctions map to increasing new slots; new suffixes are inserted at sorted SA points, some past the last entry
static Case makeCase(Kind kind) {
    Case c;
    c.kind=kind;
    c.GstrandBit=(char) (32+rng()%5);
    c.nSA1=1+rng()%(rng()%4==0 ? 20 : 3000);
    c.sjdbLength= rng()%3==0 ? 151 : 2+rng()%60;
    c.sjGstart=2+rng()%1000000;
    c.nSJ1=rng()%4==0 ? 0 : rng()%40;
    c.nSJ=c.nSJ1 + (kind==UNCHANGED ? 0 : (kind==NEW_JUNCTIONS_ONLY ? 1 : 0)+rng()%12);
    std::vector<uint32> slots(c.nSJ);
    for (uint64 i=0;i<c.nSJ;i++) slots[i]=(uint32) i;
    std::shuffle(slots.begin(), slots.end(), rng);
    c.oldSJind.assign(slots.begin(), slots.begin()+c.nSJ1);
    std::sort(c.oldSJind.begin(), c.oldSJind.end());
    if (kind==REORDERED) {//no new junctions, but two old ones swap slots
        if (c.nSJ1<2) c.nSJ=c.nSJ1=2+rng()%38;
        c.oldSJind.resize(c.nSJ1);
        for (uint64 i=0;i<c.nSJ1;i++) c.oldSJind[i]=(uint32) i;
        uint64 a=rng()%c.nSJ1, b=(a+1+rng()%(c.nSJ1-1))%c.nSJ1;
        std::swap(c.oldSJind[a], c.oldSJind[b]);
    };
    c.sjNew=c.nSJ-c.nSJ1;
    uint64 drop=0;
    if (kind==DROPPED) {//the last old junctions left out of the new set, as many new ones (if any) or none: the genome length
                        //or the number of new junctions changes while the old ones keep their slots
        if (c.nSJ1<1) c.nSJ1=1+rng()%39;
        drop=1+rng()%c.nSJ1;
        c.sjNew= rng()%2 ? drop : 0;
        c.nSJ=c.nSJ1-drop+c.sjNew;
        c.oldSJind.resize(c.nSJ1);
        for (uint64 i=0;i<c.nSJ1;i++) c.oldSJind[i]=(uint32) i;
    };

    uint64 nGenome1=c.sjGstart+c.nSJ1*c.sjdbLength, nGsj=c.nSJ*c.sjdbLength, N2bit=1LLU<<c.GstrandBit;
    for (uint64 i=0;i<c.nSA1;i++) {//no suffixes in dropped junctions, upstream cannot shift them
        bool inSj= c.nSJ1>drop && rng()%3==0;
        uint64 pos= inSj ? c.sjGstart+rng()%((c.nSJ1-drop)*c.sjdbLength) : rng()%c.sjGstart;//an old genome position
        if (rng()%2) {
            c.sa.push_back(pos);
        } else {//minus strand: the suffix of the reverse complement that ends at pos, never the one past the genome end
            if (pos==0) pos=1;
            c.sa.push_back((nGenome1-pos) | N2bit);
        };
    };

    c.nInd= kind==GENERAL && nGsj>0 ? rng()%60 : 0;
    std::vector<uint64> points(c.nInd);
    for (uint64 &p: points) p= rng()%5==0 ? c.nSA1 : rng()%(c.nSA1+1);
    std::sort(points.begin(), points.end());
    for (uint64 i=0;i<c.nInd;i++) {
        c.indArray.push_back(points[i]);
        c.indArray.push_back(rng()%(2*nGsj));
    };
    c.indArray.push_back(-999);
    c.indArray.push_back(-999);

    //spare entries past SA2's end: none or a few (SA lies right after the entries SA2 gains), room for all junctions as
    //genomeLoad and genomeGenerate leave; or room for a few 64-entry rounds, or for all entries at once
    uint64 extras[]={0, 1, 2, rng()%200, 2*c.nSJ*c.sjdbLength+rng()%10, c.nInd+64*(1+rng()%8)+rng()%64, c.nSA1+c.nInd+64+rng()%100};
    c.extra=extras[rng()%7];
    return c;
}

//the rounds must start at 0, end at the last entry, split at multiples of 64, and write each round's output before the
//first byte of SA the round reads
static bool roundsValid(const std::vector<uint64> &rounds, uint64 gapBytes, uint64 w, uint64 nSA1, const uint64 *indArray,
                        uint64 nInd) {
    if (rounds.empty()) return true;
    if (rounds.size()<2 || rounds.front()!=0 || rounds.back()!=nSA1+nInd) return false;
    for (uint64 r=0; r+1<rounds.size(); r++) {
        uint64 J0=rounds[r], J1=rounds[r+1];
        if (J1<=J0 || (r>0 && J0%64!=0)) return false;
        uint64 isa=0, isj=0;//the input of output entry J0
        for (uint64 J=0; J<J0; J++) {
            if (isj<nInd && std::min(indArray[2*isj], nSA1)==isa) isj++; else isa++;
        };
        if ((J1*w+7)/8 > gapBytes+isa*w/8) return false;
    };
    return true;
}

static bool runCase(Case &c) {
    uint64 w=c.GstrandBit+1, nSA2=c.nSA1+c.nInd;
    PackedArray buf;
    buf.defineBits(w, nSA2+c.extra);
    const uint64 guard=64;
    std::vector<char> init(buf.lengthByte+guard);
    for (char &b: init) b=(char) rng();
    PackedArray SA;
    SA.defineBits(w, c.nSA1);
    SA.pointArray(init.data()+buf.lengthByte-SA.lengthByte);
    for (uint64 i=0;i<c.nSA1;i++) SA.writePacked(i, c.sa[i]);
    std::vector<char> want(init), got(init);

    GenomeLite g, g1;
    g.nChrReal=g1.nChrReal=1;
    g.chrStart[1]=g1.chrStart[1]=c.sjGstart;
    g.sjdbLength=g1.sjdbLength=c.sjdbLength;
    g.GstrandBit=g1.GstrandBit=c.GstrandBit;
    g1.nSA=c.nSA1;
    g.nSA=nSA2;
    g1.nGenome=c.sjGstart+c.nSJ1*c.sjdbLength;
    g.nGenome=c.sjGstart+c.nSJ*c.sjdbLength;
    uint64 nGsj=c.nSJ*c.sjdbLength, sjNew=c.sjNew;

    PackedArray refSA, refSA2, candSA, candSA2;
    refSA.defineBits(w, c.nSA1); refSA2.defineBits(w, nSA2);
    candSA.defineBits(w, c.nSA1); candSA2.defineBits(w, nSA2);
    refSA.pointArray(want.data()+buf.lengthByte-SA.lengthByte); refSA2.pointArray(want.data());
    candSA.pointArray(got.data()+buf.lengthByte-SA.lengthByte); candSA2.pointArray(got.data());
    std::vector<uint32> oldSJindRef(c.oldSJind);
    oldSJindRef.push_back(0);//keeps .data() valid when there are no old junctions
    std::vector<uint64> indArrayRef(c.indArray), indArrayCand(c.indArray);
    ref::insertSA(refSA, refSA2, g, g1, oldSJindRef.data(), indArrayRef.data(), c.nInd, nGsj, sjNew);
    const uint64 threadChoices[]={1, 2, 3, 4, 7};
    uint64 nThreads=threadChoices[rng()%5], gapBytes=buf.lengthByte-SA.lengthByte;
    SAinsertReport rep=sjdbInsertSAindices(candSA, candSA2, c.nSA1, g1.nGenome, g.nGenome, c.sjGstart, c.sjdbLength, c.nSJ1,
                                           oldSJindRef.data(), sjNew, nGsj, indArrayCand.data(), c.nInd, c.GstrandBit, nThreads, 1);
    bool block=rep.block;

    bool ok=true;
    for (uint64 i=0;i<want.size();i++) {
        if (want[i]!=got[i]) {
            if (nFail++<5) printf("  kind %d nSA1 %llu nInd %llu nSJ1 %llu nSJ %llu extra %llu w %llu: byte %llu of %llu differs\n",
                                  c.kind, c.nSA1, c.nInd, c.nSJ1, c.nSJ, c.extra, w, i, (uint64) want.size());
            ok=false;
            break;
        };
    };
    bool unchanged= c.nInd==0 && c.sjNew==0 && c.nSJ==c.nSJ1;
    for (uint64 i=0;i<c.nSJ1;i++) unchanged= unchanged && c.oldSJind[i]==i;
    if (c.kind==UNCHANGED && !unchanged) {
        printf("  kind UNCHANGED made a changing case\n");
        nFail++;
    };
    if (block != unchanged) {
        if (nFail++<5) printf("  kind %d nSA1 %llu nInd %llu nSJ1 %llu nSJ %llu: block move %d\n", c.kind, c.nSA1, c.nInd, c.nSJ1, c.nSJ, block);
        ok=false;
    };
    //a block move is shared when the room holds a round of nThreads entries; it then takes rounds of the room's size
    uint64 nBytes=(c.nSA1-1)*w/8+1, wantRounds=0;
    if (block) {
        if (nThreads>1 && gapBytes*8/w>=nThreads) wantRounds=(nBytes+gapBytes-1)/gapBytes;
        nBlockShared+= wantRounds>0;
    } else {
        std::vector<uint64> rounds=sjdbInsertSArounds(gapBytes, w, c.nSA1, c.indArray.data(), c.nInd, nThreads);
        if (!roundsValid(rounds, gapBytes, w, c.nSA1, c.indArray.data(), c.nInd)) {
            if (nFail++<5) printf("  kind %d nSA1 %llu nInd %llu gap %llu w %llu: invalid rounds\n", c.kind, c.nSA1, c.nInd, gapBytes, w);
            ok=false;
        };
        if (nThreads>1 && !rounds.empty()) wantRounds=rounds.size()-1;
        nMultiRound+= wantRounds>1;
    };
    if (rep.nRounds!=wantRounds) {
        if (nFail++<5) printf("  kind %d nSA1 %llu nInd %llu gap %llu threads %llu: %llu rounds, want %llu\n", c.kind, c.nSA1, c.nInd,
                              gapBytes, nThreads, rep.nRounds, wantRounds);
        ok=false;
    };
    return ok;
}

static bool testMatchesUpstream(Kind kind, int nCases) {
    nFail=0;
    for (int n=0;n<nCases;n++) {
        Case c=makeCase(kind);
        runCase(c);
    };
    return nFail==0;
}

//rounds for large arrays and many insertions, with room from barely enough to plenty
static bool testRoundsBigArrays() {
    nFail=0;
    uint64 nPlanned=0;
    for (int n=0;n<300;n++) {
        uint64 w=33+rng()%5, nSA1=1+rng()%2000000, nInd=rng()%20000;
        std::vector<uint64> ind;
        std::vector<uint64> points(nInd);
        for (uint64 &p: points) p=rng()%(nSA1+1);
        std::sort(points.begin(), points.end());
        for (uint64 p: points) {ind.push_back(p); ind.push_back(0);};
        ind.push_back(-999); ind.push_back(-999);
        uint64 gapBytes=nInd*w/8+16+rng()%(1+(rng()%2 ? 4000 : nSA1*w/8));
        uint64 minRound=1+rng()%5000;
        std::vector<uint64> rounds=sjdbInsertSArounds(gapBytes, w, nSA1, ind.data(), nInd, minRound);
        nPlanned+= rounds.size()>2;
        if (!roundsValid(rounds, gapBytes, w, nSA1, ind.data(), nInd)) {
            if (nFail++<5) printf("  nSA1 %llu nInd %llu gap %llu w %llu: invalid rounds\n", nSA1, nInd, gapBytes, w);
        };
        for (uint64 r=1; r+2<rounds.size(); r++) {//all but the last round are long enough
            if (rounds[r+1]-rounds[r]<minRound) {
                if (nFail++<5) printf("  round %llu of %llu entries, shorter than %llu\n", r, rounds[r+1]-rounds[r], minRound);
                break;
            };
        };
    };
    if (nPlanned<50) {printf("  only %llu cases planned 2+ rounds\n", nPlanned); nFail++;};
    return nFail==0;
}

int main() {
    int rc=0;
    if (!testRoundsBigArrays()) { printf("FAIL testRoundsBigArrays (%d)\n", nFail); rc=1; }
    if (!testMatchesUpstream(UNCHANGED, 3000)) { printf("FAIL testMatchesUpstream(UNCHANGED) (%d)\n", nFail); rc=1; }
    if (!testMatchesUpstream(GENERAL, 3000)) { printf("FAIL testMatchesUpstream(GENERAL) (%d)\n", nFail); rc=1; }
    if (!testMatchesUpstream(NEW_JUNCTIONS_ONLY, 3000)) { printf("FAIL testMatchesUpstream(NEW_JUNCTIONS_ONLY) (%d)\n", nFail); rc=1; }
    if (!testMatchesUpstream(REORDERED, 3000)) { printf("FAIL testMatchesUpstream(REORDERED) (%d)\n", nFail); rc=1; }
    if (!testMatchesUpstream(DROPPED, 3000)) { printf("FAIL testMatchesUpstream(DROPPED) (%d)\n", nFail); rc=1; }
    if (nMultiRound<300 || nBlockShared<100) {
        printf("FAIL coverage: %llu cases rebuilt in 2+ rounds, %llu shared block moves\n", nMultiRound, nBlockShared); rc=1;
    };
    return rc;
}
