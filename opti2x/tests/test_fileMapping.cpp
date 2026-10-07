// mapFileReadOnly must give the bytes a private copy would hold: the file's bytes from offset, then zeros where the
// requested length runs past the end of the file (PackedArray::allocateArray zeroes the last word, which the SA and
// SAindex files can be a few bytes short of). Covers file sizes at, below and above page boundaries.
#include "fileMapping.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

static std::string dir;
static int nFail=0;

static std::string writeFile(const std::vector<char>& bytes) {
    static int n=0;
    std::string path=dir+"/f"+std::to_string(n++);
    FILE* f=fopen(path.c_str(), "wb");
    if (!f || fwrite(bytes.data(), 1, bytes.size(), f)!=bytes.size()) { printf("cannot write %s\n", path.c_str()); exit(2); }
    fclose(f);
    return path;
}

static bool testMappedBytesMatchFileThenZeros() {
    nFail=0;
    std::mt19937_64 rng(20261007);
    uint64 pg=sysconf(_SC_PAGESIZE);
    for (uint64 size: {(uint64) 1, (uint64) 100, pg-1, pg, pg+1, 3*pg, 3*pg+17, 1000003LLU}) {
        std::vector<char> bytes(size);
        for (char& c: bytes) c=(char) (1+rng()%255);
        std::string path=writeFile(bytes);
        for (uint64 offset: {(uint64) 0, (uint64) 8, (uint64) 1000, pg+3}) {
            if (offset>=size) continue;
            for (uint64 extra: {(uint64) 0, (uint64) 1, (uint64) 8, pg, 2*pg+5}) {
                uint64 inFile=size-offset, length=inFile+extra;
                for (uint64 len: {length, inFile/2+1}) {
                    char* m=mapFileReadOnly(path, offset, len);
                    if (m==NULL) { if (nFail++<5) printf("  size %llu offset %llu length %llu: NULL\n", size, offset, len); continue; }
                    for (uint64 i=0;i<len;i++) {
                        char want= i<inFile ? bytes[offset+i] : 0;
                        if (m[i]!=want) {
                            if (nFail++<5) printf("  size %llu offset %llu length %llu: byte %llu is %d, want %d\n",
                                                  size, offset, len, i, m[i], want);
                            break;
                        }
                    }
                }
            }
        }
    }
    return nFail==0;
}

static bool testMissingFileGivesNull() {
    nFail=0;
    if (mapFileReadOnly(dir+"/missing", 0, 100)!=NULL) { printf("  mapped a missing file\n"); nFail++; }
    return nFail==0;
}

int main() {
    const char* tmp=getenv("TMPDIR");
    std::string base=std::string(tmp ? tmp : "/tmp")+"/star_fileMapping_XXXXXX";
    std::vector<char> buf(base.begin(), base.end()); buf.push_back(0);
    if (!mkdtemp(buf.data())) { printf("mkdtemp failed\n"); return 2; }
    dir=buf.data();
    int rc=0;
    if (!testMappedBytesMatchFileThenZeros()) { printf("FAIL testMappedBytesMatchFileThenZeros (%d)\n", nFail); rc=1; }
    if (!testMissingFileGivesNull()) { printf("FAIL testMissingFileGivesNull (%d)\n", nFail); rc=1; }
    std::string rm="rm -rf "+dir;
    if (system(rm.c_str())!=0) printf("could not remove %s\n", dir.c_str());
    return rc;
}
