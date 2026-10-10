// Differential test: the candidate stitchAlignToTranscript must return the same score and leave the transcript in the
// same state as upstream STAR 2.7.11b on randomized reads stitched across junctions, deletions, insertions and mates.
// Upstream stitchAlignToTranscript and the extendAlign it calls are compiled into namespace ref (UPSTREAM_STITCH and
// UPSTREAM_EXTENDALIGN, the .cpp files); extendAlign is renamed there, since argument-dependent lookup on Transcript*
// would otherwise make the call ambiguous with the candidate's ::extendAlign.
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Genome.h"
#include "Transcript.h"
#include "extendAlign.h"
#include "binarySearch2.h"
#include "stitchAlignToTranscript.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <map>
#include <random>

namespace ref {
#define extendAlign upstreamExtendAlign
#include UPSTREAM_EXTENDALIGN
#include UPSTREAM_STITCH
#undef extendAlign
}

static std::mt19937_64 rng(20261007);
static uint rnd(uint n) {return rng()%n;};             //0..n-1
static bool chance(double p) {return (rng()%1000000)<p*1000000;};

static const uint gMargin=1000, gLen=8000, rMargin=700, rLen=600;

struct Case {
    vector<char> gBuf, rBuf;
    vector<uint> sjS, sjE; vector<uint8> sjMotif, sjShL, sjShR, sjStr;
    uint rAend, gAend, rBstart, gBstart, L, iFragB, sjAB, mmMaxTotal;
    Transcript tr;
};

static const char motifs[6][4]={{2,3,0,2},{1,3,0,1},{2,1,0,2},{1,3,2,1},{0,3,0,1},{2,3,0,3}};

static void setParameters(Parameters &P) {
    P.alignIntronMin = chance(0.8) ? 21 : 5+rnd(60);
    P.alignIntronMax = chance(0.7) ? 0 : 50+rnd(3000);
    P.alignMatesGapMax = chance(0.7) ? 0 : 100+rnd(2000);
    P.alignSJstitchMismatchNmax = {0,-1,0,0};
    if (chance(0.3)) for (auto &x: P.alignSJstitchMismatchNmax) x=(int32)rnd(4)-1;
    P.alignEndsProtrude.nBasesMax = chance(0.8) ? 0 : rnd(20);
    P.alignInsertionFlush.flushRight = chance(0.5);
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) P.alignEndsType.ext[i][j] = chance(0.15);
    P.outFilterMismatchNoverLmax = chance(0.7) ? 0.3 : 0.02+0.98*rnd(1000)/1000.0;
    P.scoreGap = chance(0.8) ? 0 : -(intScore)rnd(4);
    //motif penalties may be set positive: the junction scan's early stop must bound them too
    P.scoreGapNoncan = chance(0.8) ? -8 : (intScore)rnd(16)-12;
    P.scoreGapGCAG = chance(0.8) ? -4 : (intScore)rnd(12)-8;
    P.scoreGapATAC = chance(0.8) ? -8 : (intScore)rnd(16)-12;
    P.scoreDelBase = -2; P.scoreDelOpen = -2; P.scoreInsBase = -2; P.scoreInsOpen = -2;
    if (chance(0.2)) {P.scoreDelBase=-(intScore)rnd(4); P.scoreDelOpen=-(intScore)rnd(4); P.scoreInsBase=-(intScore)rnd(4); P.scoreInsOpen=-(intScore)rnd(4);};
    P.scoreStitchSJshift = chance(0.8) ? 1 : rnd(4);
    P.pGe.sjdbScore = chance(0.8) ? 2 : rnd(5);
};

static void mutate(char *s, uint n, double mmRate) {
    for (uint i=0;i<n;i++) {
        if (chance(mmRate)) s[i]=(s[i]+1+rnd(3))%4;
        if (chance(0.003)) s[i]=4;
    };
};

//builds a genome with junctions/indels planted between "exon A" and "exon B", the read spanning them, seeds A (ending at
//rAend/gAend, already in the transcript) and B (rBstart/gBstart/L), and annotated junctions
static void makeCase(Case &c, Parameters &P) {
    c.gBuf.assign(gLen+2*gMargin, GENOME_spacingChar);
    char *G=c.gBuf.data()+gMargin;
    for (uint i=0;i<gLen;i++) G[i]= chance(0.002) ? 4 : rnd(4);
    if (chance(0.1)) {uint s=rnd(gLen-200); for (uint i=0;i<100;i++) G[s+i]=G[s+(i%(1+rnd(4)))];}; //short tandem repeats

    c.rBuf.assign(rLen+2*rMargin, MARK_FRAG_SPACER_BASE);
    char *R=c.rBuf.data()+rMargin;

    const int kind=rnd(10); //0-4 junction/deletion, 5-6 insertion, 7 no gap, 8 mates, 9 sjdb fast path
    uint LA=12+rnd(80), LB=12+rnd(80);
    uint gA0=300+rnd(1500), rA0=rnd(60);
    uint gapG=0, gapR=0; //bases between exon A and exon B in genome / read
    if (kind<=4 || kind==9) {
        gapG = chance(0.3) ? 1+rnd(P.alignIntronMin+3) : 1+rnd(3000);
    } else if (kind<=6) {
        gapR=1+rnd(8);
        if (chance(0.3)) gapG=1+rnd(gapR+3);
    };
    uint gB0=gA0+LA+gapG;
    if (gB0+LB+300>gLen) {gB0=gLen-LB-300;};

    if (gapG>=4 && chance(0.7)) {//plant a motif at the intron ends
        const char *m=motifs[rnd(6)];
        G[gA0+LA]=m[0]; G[gA0+LA+1]=m[1]; G[gB0-2]=m[2]; G[gB0-1]=m[3];
    };

    uint rB0=rA0+LA+gapR;
    memcpy(R+rA0, G+gA0, LA);
    for (uint i=0;i<gapR;i++) R[rA0+LA+i]=rnd(4);
    memcpy(R+rB0, G+gB0, LB);
    uint Lread=rB0+LB+rnd(30);
    for (uint i=rB0+LB;i<Lread;i++) R[i]=rnd(4);
    for (uint i=0;i<rA0;i++) R[i]=rnd(4);
    mutate(R, Lread, chance(0.5) ? 0.0 : 0.06*rnd(1000)/1000.0);

    //seed A: last exon of the transcript, ends near the true exon A end
    int sA = (int)rnd(14)-4;
    c.rAend = rA0+LA-1-sA; if (c.rAend<rA0+2) c.rAend=rA0+2;
    c.gAend = gA0+(c.rAend-rA0);
    //seed B: starts near the true exon B start
    int sB = (int)rnd(16)-6;
    c.rBstart = rB0+sB; if (c.rBstart<rA0+1) c.rBstart=rA0+1;
    c.gBstart = gB0+(int)(c.rBstart-rB0);
    uint rBmax=rB0+LB-1;
    c.L = (c.rBstart<=rBmax) ? 1+rnd(rBmax-c.rBstart+1) : 1+rnd(5);
    c.iFragB=0;
    c.sjAB=(uint)-1;

    uint iFragA=0;
    if (kind==7) {//contiguous seeds on the same diagonal
        c.rBstart=c.rAend+1+rnd(10); c.gBstart=c.gAend+(c.rBstart-c.rAend);
    } else if (kind==8) {//mates: B on fragment 1, downstream in the genome
        uint spacer=c.rAend+1+rnd(30);
        if (spacer>=rLen-120) spacer=rLen-120;
        for (uint i=c.rAend+1;i<spacer;i++) R[i]=G[c.gAend+1+(i-c.rAend-1)];
        R[spacer]=MARK_FRAG_SPACER_BASE;
        uint m2=spacer+1, Lm2=40+rnd(70);
        uint gM2=c.gAend+(chance(0.2) ? -(int)rnd(80) : (int)rnd(600));
        memcpy(R+m2, G+gM2, Lm2);
        mutate(R+m2, Lm2, 0.02);
        c.rBstart=m2+rnd(20); c.gBstart=gM2+(c.rBstart-m2); c.L=1+rnd(m2+Lm2-c.rBstart);
        c.iFragB=1;
    };

    //annotated junctions: the planted junction (sometimes) and random decoys, sorted by start then end
    struct Sj {uint s, e; uint8 m, l, r, str;};
    vector<Sj> sj;
    uint jS=gA0+LA, jE=gB0-1;
    if (gapG>0 && chance(0.5)) sj.push_back({jS, jE, (uint8)rnd(7), (uint8)(chance(0.5)?0:rnd(4)), (uint8)(chance(0.5)?0:rnd(4)), (uint8)rnd(3)});
    if (gapG>0 && chance(0.3)) sj.push_back({jS+rnd(5), jE+rnd(5), (uint8)rnd(7), 0, 0, (uint8)rnd(3)}); //near-miss decoy
    for (uint n=rnd(6), i=0; i<n; i++) {uint s=rnd(gLen-500); sj.push_back({s, s+20+rnd(400), (uint8)rnd(7), (uint8)rnd(3), (uint8)rnd(3), (uint8)rnd(3)});};
    std::sort(sj.begin(), sj.end(), [](const Sj &a, const Sj &b){return a.s<b.s || (a.s==b.s && a.e<b.e);});
    c.sjS.clear(); c.sjE.clear(); c.sjMotif.clear(); c.sjShL.clear(); c.sjShR.clear(); c.sjStr.clear();
    for (auto &x: sj) {c.sjS.push_back(x.s); c.sjE.push_back(x.e); c.sjMotif.push_back(x.m); c.sjShL.push_back(x.l); c.sjShR.push_back(x.r); c.sjStr.push_back(x.str);};

    //transcript: optional earlier exons, last exon = seed A
    c.tr.reset();
    uint nEx=1+(chance(0.3) ? rnd(3) : 0);
    uint exL=c.rAend-rA0+1; if (chance(0.3)) exL=1+rnd(exL);
    uint exR=c.rAend-exL+1;
    for (uint i=0;i<nEx;i++) {
        uint iex=nEx-1-i;
        c.tr.exons[iex][EX_R]= (i==0) ? exR : 0;
        c.tr.exons[iex][EX_G]= (i==0) ? c.gAend-exL+1 : gA0-400+rnd(100);
        c.tr.exons[iex][EX_L]= (i==0) ? exL : 5+rnd(20);
        c.tr.exons[iex][EX_iFrag]=iFragA;
        c.tr.exons[iex][EX_sjA]=(uint)-1;
        if (iex<nEx-1) {c.tr.canonSJ[iex]=(int)rnd(8)-1; c.tr.sjAnnot[iex]=rnd(2); c.tr.sjStr[iex]=rnd(3); c.tr.shiftSJ[iex][0]=rnd(3); c.tr.shiftSJ[iex][1]=rnd(3);};
    };
    c.tr.nExons=nEx;
    c.tr.rStart=c.tr.exons[0][EX_R]; c.tr.gStart=c.tr.exons[0][EX_G];
    c.tr.nMM=rnd(4); c.tr.nMatch=exL+rnd(30); c.tr.maxScore=rnd(100);
    if (chance(0.2)) {c.tr.nGap=rnd(2); c.tr.lGap=rnd(500); c.tr.nDel=rnd(2); c.tr.lDel=rnd(5); c.tr.nIns=rnd(2); c.tr.lIns=rnd(5);};

    if (kind==9 && !c.sjS.empty()) {//annotated junction carried by seed A and B adjacent in the read
        uint i=rnd(c.sjS.size());
        c.sjAB=i;
        c.tr.exons[nEx-1][EX_sjA]=i;
        c.rBstart=c.rAend+1;
        c.gBstart= chance(0.8) ? c.sjE[i]+1+(c.gAend+1>c.sjS[i] ? 0 : 0) : c.gAend+1+rnd(3);
        if (c.gBstart+c.L+10>gLen) c.gBstart=gLen-c.L-10;
    } else if (chance(0.05)) {
        c.sjAB=rnd(4); //sjA that does not match the transcript
    };

    c.mmMaxTotal = chance(0.8) ? 10 : rnd(12);
    if (chance(0.03)) c.tr.nExons=MAX_N_EXONS; //exon limit guard
};

static std::map<int,uint> outcomes;
static uint nFail=0;

static bool caseFailed;
static void report(const char *what, uint iCase, long long u, long long c) {
    if (nFail<8) printf("  case %llu: %s upstream %lld candidate %lld\n", iCase, what, u, c);
    caseFailed=true;
};

static void compare(uint iCase, intScore su, intScore sc, const Transcript &u, const Transcript &c) {
    caseFailed=false;
    if (su!=sc) report("score", iCase, su, sc);
    #define CMP(x) if (u.x!=c.x) report(#x, iCase, (long long)u.x, (long long)c.x);
    CMP(nExons) CMP(nMM) CMP(nMatch) CMP(nGap) CMP(lGap) CMP(nDel) CMP(lDel) CMP(nIns) CMP(lIns) CMP(maxScore) CMP(nUnique) CMP(extendL)
    #undef CMP
    if (memcmp(u.exons, c.exons, sizeof(u.exons))) report("exons differ", iCase, 0, 1);
    if (memcmp(u.shiftSJ, c.shiftSJ, sizeof(u.shiftSJ))) report("shiftSJ differ", iCase, 0, 1);
    if (memcmp(u.canonSJ, c.canonSJ, sizeof(u.canonSJ))) report("canonSJ differ", iCase, 0, 1);
    if (memcmp(u.sjAnnot, c.sjAnnot, sizeof(u.sjAnnot))) report("sjAnnot differ", iCase, 0, 1);
    if (memcmp(u.sjStr, c.sjStr, sizeof(u.sjStr))) report("sjStr differ", iCase, 0, 1);
    if (caseFailed) nFail++;
};

int main(int argc, char **argv) {
    uint nCases = argc>1 ? atoi(argv[1]) : 1000000;
    Parameters P;
    Genome mapGen(P, P.pGe);
    Case c;
    for (uint iCase=0; iCase<nCases; iCase++) {
        setParameters(P);
        makeCase(c, P);
        mapGen.G=c.gBuf.data()+gMargin;
        mapGen.nGenome=gLen;
        mapGen.sjdbN=c.sjS.size();
        mapGen.sjdbStart=c.sjS.data(); mapGen.sjdbEnd=c.sjE.data(); mapGen.sjdbMotif=c.sjMotif.data();
        mapGen.sjdbShiftLeft=c.sjShL.data(); mapGen.sjdbShiftRight=c.sjShR.data(); mapGen.sjdbStrand=c.sjStr.data();
        if (chance(0.2)) mapGen.sjdbN=0;
        char *R=c.rBuf.data()+rMargin;

        Transcript tu, tc;
        for (Transcript *t: {&tu, &tc}) {//identical contents beyond nExons too, so whole arrays can be compared
            memset(t->exons, 0, sizeof(t->exons)); memset(t->shiftSJ, 0, sizeof(t->shiftSJ)); memset(t->canonSJ, 0, sizeof(t->canonSJ));
            memset(t->sjAnnot, 0, sizeof(t->sjAnnot)); memset(t->sjStr, 0, sizeof(t->sjStr));
        };
        tu.copyStitchCore(c.tr); tu.maxScore=c.tr.maxScore; tu.nUnique=c.tr.nUnique; tu.extendL=c.tr.extendL;
        tc.copyStitchCore(c.tr); tc.maxScore=c.tr.maxScore; tc.nUnique=c.tr.nUnique; tc.extendL=c.tr.extendL;
        intScore su=ref::stitchAlignToTranscript(c.rAend, c.gAend, c.rBstart, c.gBstart, c.L, c.iFragB, c.sjAB, P, R, mapGen, &tu, c.mmMaxTotal);
        intScore sc=stitchAlignToTranscript(c.rAend, c.gAend, c.rBstart, c.gBstart, c.L, c.iFragB, c.sjAB, P, R, mapGen, &tc, c.mmMaxTotal);
        compare(iCase, su, sc, tu, tc);
        int key = su<-1000000 ? su : (tu.nExons>c.tr.nExons ? 100+tu.canonSJ[tu.nExons-2] : 0);
        outcomes[key]++;
    };
    fprintf(stderr, "outcomes (return code, or 100+canonSJ of the new junction, or 0 = merged into the last exon):\n");
    for (auto &o: outcomes) fprintf(stderr, "  %d: %llu\n", o.first, o.second);
    if (nFail>0) {printf("FAIL: %llu of %llu cases differ\n", nFail, nCases); return 1;};
    return 0;
};
