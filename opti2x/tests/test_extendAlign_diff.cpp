// Differential fuzz test: extendAlign in source/ must return exactly what upstream STAR 2.7.11b's
// extendAlign returns. The upstream file is compiled into namespace ref (its own includes are
// no-ops here because the headers are already included).
//
// Build/run: opti2x/tests/run_unit_tests.sh
#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "extendAlign.h"

namespace ref {
#include UPSTREAM_EXTENDALIGN
}

#include <cstdio>
#include <random>

static const int PAD = 512;

int main() {
    std::mt19937_64 rng(20261006);
    const double pMMmaxChoices[] = {0.04, 0.05, 0.066, 0.1, 0.13, 0.3, 1.0};
    char rBuf[2*PAD], gBuf[2*PAD];
    uint64 nFail = 0, nCases = 2000000;

    for (uint64 iCase = 0; iCase < nCases; iCase++) {
        // mostly matches with a tunable mismatch rate, occasional N (4), chr padding (5) and spacer
        double pMM = std::uniform_real_distribution<double>(0.0, 0.4)(rng);
        for (int i = 0; i < 2*PAD; i++) {
            char b = rng() % 4;
            rBuf[i] = b;
            gBuf[i] = (std::uniform_real_distribution<double>(0, 1)(rng) < pMM) ? char((b + 1 + rng() % 3) % 4) : b;
            uint r = rng() % 1000;
            if (r < 4) rBuf[i] = 4;
            else if (r < 8) gBuf[i] = 4;
            else if (r < 9) gBuf[i] = 5;
            else if (r < 10) rBuf[i] = MARK_FRAG_SPACER_BASE;
        };

        int dir = (rng() % 2) ? 1 : -1;
        uint L = 1 + rng() % 200;
        uint Lprev = rng() % 150;
        uint nMMprev = rng() % 6;
        uint nMMmax = rng() % 12;
        double pMMmax = pMMmaxChoices[rng() % (sizeof(pMMmaxChoices)/sizeof(pMMmaxChoices[0]))];
        bool extendToEnd = (rng() % 8) == 0;
        // gStart near 0 exercises the (gStart+iG)==(uint)-1 guard when extending left
        uint gStart = (rng() % 16 == 0) ? rng() % 8 : 256;
        char *G = gBuf + PAD - gStart;
        char *R = rBuf;
        uint rStart = PAD;

        Transcript trCand, trRef;
        trCand.extendL = trRef.extendL = 7777;
        trCand.maxScore = trRef.maxScore = 7777;
        trCand.nMatch = trRef.nMatch = 7777;
        trCand.nMM = trRef.nMM = 7777;

        bool yCand = extendAlign(R, G, rStart, gStart, dir, dir, L, Lprev, nMMprev, nMMmax, pMMmax, extendToEnd, &trCand);
        bool yRef = ref::extendAlign(R, G, rStart, gStart, dir, dir, L, Lprev, nMMprev, nMMmax, pMMmax, extendToEnd, &trRef);

        if (yCand != yRef || trCand.extendL != trRef.extendL || trCand.maxScore != trRef.maxScore
                || trCand.nMatch != trRef.nMatch || trCand.nMM != trRef.nMM) {
            if (nFail < 5)
                printf("MISMATCH case %llu: dir=%d L=%llu Lprev=%llu nMMprev=%llu nMMmax=%llu pMMmax=%g end=%d | "
                       "cand ret=%d extL=%llu score=%d nMM=%llu | ref ret=%d extL=%llu score=%d nMM=%llu\n",
                       (unsigned long long) iCase, dir, (unsigned long long) L, (unsigned long long) Lprev,
                       (unsigned long long) nMMprev, (unsigned long long) nMMmax, pMMmax, extendToEnd,
                       yCand, (unsigned long long) trCand.extendL, trCand.maxScore, (unsigned long long) trCand.nMM,
                       yRef, (unsigned long long) trRef.extendL, trRef.maxScore, (unsigned long long) trRef.nMM);
            nFail++;
        };
    };

    printf("test_extendAlign_diff: %llu cases, %llu mismatches\n", (unsigned long long) nCases, (unsigned long long) nFail);
    return nFail == 0 ? 0 : 1;
}
