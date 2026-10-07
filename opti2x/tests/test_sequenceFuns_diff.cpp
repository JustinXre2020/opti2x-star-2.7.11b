// Differential test: the per-read sequence conversions must produce exactly the bytes upstream STAR 2.7.11b produces,
// for every input byte value. Upstream SequenceFuns.cpp (UPSTREAM_SEQUENCEFUNS) is compiled into namespace ref.
#include "SequenceFuns.h"
#include <assert.h>
#include <cstdio>
#include <cstring>
#include <random>

namespace ref {
#include UPSTREAM_SEQUENCEFUNS
}

static std::mt19937_64 rng(20261007);
static int nFail=0;

static void fail(const char *test, const vector<char> &in, uint L) {
    if (nFail++<5) {
        printf("  %s differs, L=%llu, input bytes:", test, (unsigned long long)L);
        for (uint i=0; i<L && i<16; i++) printf(" %d", in[i]);
        printf("\n");
    };
}

//inputs: every byte value at every offset of a short read, then random reads, mostly sequence letters
static vector<vector<char>> inputs() {
    vector<vector<char>> v;
    for (int b=0; b<256; b++) for (uint L: {1, 7, 16, 33}) {
        vector<char> s(L);
        for (uint i=0; i<L; i++) s[i]=(char)((b+i*37)%256);
        v.push_back(s);
    };
    const char letters[]="ACGTNacgtnRYKMSWBDVHrykmswbdvh=.-\0\1\2\3\4\5";
    for (int n=0; n<20000; n++) {
        vector<char> s(rng()%700);
        bool raw=rng()%4==0;
        for (char &c: s) c= raw ? (char)rng() : letters[rng()%(sizeof(letters)-1)];
        v.push_back(s);
    };
    v.push_back(vector<char>());
    return v;
}

static bool testComplementSeqNumbers() {
    nFail=0;
    for (auto &in: inputs()) {
        uint L=in.size();
        vector<char> u(L+8,'#'), c(L+8,'#');
        ref::complementSeqNumbers((char*)in.data(), u.data(), L);
        complementSeqNumbers((char*)in.data(), c.data(), L);
        if (u!=c) fail("complementSeqNumbers", in, L);
    };
    return nFail==0;
}

static bool testConvertNucleotidesToNumbers() {
    nFail=0;
    for (auto &in: inputs()) {
        uint L=in.size();
        vector<char> u(L+8,'#'), c(L+8,'#');
        ref::convertNucleotidesToNumbers(in.data(), u.data(), L);
        convertNucleotidesToNumbers(in.data(), c.data(), L);
        if (u!=c) fail("convertNucleotidesToNumbers", in, L);
    };
    return nFail==0;
}

static bool testRevComplementNucleotides() {
    nFail=0;
    for (auto &in: inputs()) {
        uint L=in.size();
        vector<char> u(L+8,'#'), c(L+8,'#');
        ref::revComplementNucleotides((char*)in.data(), u.data(), L);
        revComplementNucleotides((char*)in.data(), c.data(), L);
        if (u!=c) fail("revComplementNucleotides(char*)", in, L);
        string su(in.begin(), in.end()), sc(su);
        ref::revComplementNucleotides(su);
        revComplementNucleotides(sc);
        if (su!=sc) fail("revComplementNucleotides(string&)", in, L);
    };
    return nFail==0;
}

int main() {
    int rc=0;
    if (!testComplementSeqNumbers()) {printf("FAIL testComplementSeqNumbers (%d)\n", nFail); rc=1;};
    if (!testConvertNucleotidesToNumbers()) {printf("FAIL testConvertNucleotidesToNumbers (%d)\n", nFail); rc=1;};
    if (!testRevComplementNucleotides()) {printf("FAIL testRevComplementNucleotides (%d)\n", nFail); rc=1;};
    return rc;
}
