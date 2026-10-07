// Unit test for WinBinReset: per-read reset of the winBin table that touches only the entries written
// since the previous reset, instead of memsetting 2*winBinN entries for every read.
//
// Build/run: opti2x/tests/run_unit_tests.sh test_winBinReset
#include "IncludeDefine.h"
#include "WinBinReset.h"

#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

static int nFail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { nFail++; printf("FAIL line %d: ", __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

struct Table {
    uint N;
    std::vector<uintWinBin> s0, s1;
    uintWinBin *wb[2];
    explicit Table(uint n) : N(n), s0(n, uintWinBinMax), s1(n, uintWinBinMax) { wb[0] = s0.data(); wb[1] = s1.data(); }
    // write value v to bins [b1,b2] of strand str, and tell the tracker about it, as STAR's window code does
    void write(WinBinReset &r, uint str, uint b1, uint b2, uintWinBin v) {
        for (uint b = b1; b <= b2; b++) wb[str][b] = v;
        r.mark(str, b1, b2);
    }
    uint64 nonEmpty() const {
        uint64 n = 0;
        for (uint i = 0; i < N; i++) n += (s0[i] != uintWinBinMax) + (s1[i] != uintWinBinMax);
        return n;
    }
};

static void testResetRestoresAllMarkedEntries() {
    std::mt19937_64 rng(1);
    for (int iter = 0; iter < 20000; iter++) {
        Table t(5000 + rng() % 5000);
        WinBinReset r;
        r.init(t.wb, t.N);
        for (int read = 0; read < 4; read++) {
            uint nWrites = rng() % 40;
            for (uint i = 0; i < nWrites; i++) {
                uint str = rng() % 2;
                uint b1 = rng() % t.N;
                uint b2 = std::min<uint>(t.N - 1, b1 + rng() % 30);
                t.write(r, str, b1, b2, rng() % 1000);
            }
            r.reset();
            CHECK(t.nonEmpty() == 0, "iter %d read %d: %llu entries left set after reset", iter, read, (unsigned long long) t.nonEmpty());
        }
    }
}

static void testResetLeavesUnmarkedEntriesAlone() {
    // a sentinel written without mark() survives a reset that did not overflow: the reset is not a full memset
    Table t(100000);
    WinBinReset r;
    r.init(t.wb, t.N);
    t.write(r, 0, 10, 20, 3);
    t.write(r, 1, 90000, 90005, 4);
    t.wb[0][50000] = 7;
    r.reset();
    CHECK(t.wb[0][15] == uintWinBinMax && t.wb[1][90003] == uintWinBinMax, "marked entries not reset");
    CHECK(t.wb[0][50000] == 7, "unmarked entry was reset (full memset instead of partial)");
}

static void testAdjacentMarksCoalesce() {
    // windows grow bin by bin; touching intervals on one strand must not use up the interval list
    Table t(1000000);
    WinBinReset r;
    r.init(t.wb, t.N);
    for (uint b = 1000; b < 1000 + 10 * WinBinReset::nMax; b++)
        t.write(r, 0, b, b, 1);
    t.wb[0][500000] = 7;
    r.reset();
    CHECK(t.nonEmpty() == 1 && t.wb[0][500000] == 7, "adjacent marks overflowed the list or were not reset");
}

static void testOverflowFallsBackToFullReset() {
    Table t(1000000);
    WinBinReset r;
    r.init(t.wb, t.N);
    for (uint i = 0; i < WinBinReset::nMax + 5; i++)
        t.write(r, i % 2, 100 + i * 50, 100 + i * 50 + 3, 2); // disjoint, alternating strands
    t.wb[0][999000] = 7; // unmarked: only a full reset clears it
    r.reset();
    CHECK(t.nonEmpty() == 0, "overflowed list did not fall back to a full reset");
    // and the tracker is usable again afterwards, in partial mode
    t.write(r, 1, 5, 6, 1);
    t.wb[1][999000] = 7;
    r.reset();
    CHECK(t.wb[1][5] == uintWinBinMax && t.wb[1][999000] == 7, "tracker not back in partial mode after an overflow");
}

int main() {
    testResetRestoresAllMarkedEntries();
    testResetLeavesUnmarkedEntriesAlone();
    testAdjacentMarksCoalesce();
    testOverflowFallsBackToFullReset();
    printf("test_winBinReset: %d failures\n", nFail);
    return nFail == 0 ? 0 : 1;
}
