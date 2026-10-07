#include "fileMapping.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

char* mapFileReadOnly(const string& path, uint64 offset, uint64 length) {
    int fd=open(path.c_str(), O_RDONLY);
    if (fd<0)
        return NULL;
    struct stat st;
    if (fstat(fd, &st)!=0) {
        close(fd);
        return NULL;
    };
    uint64 end=offset+length, fileLen=min((uint64) st.st_size, end);
    //reserve zero pages for the whole range, then lay the file over its start: pages past the end of the file stay zero
    char* base=(char*) mmap(NULL, end, PROT_READ, MAP_PRIVATE|MAP_ANON, -1, 0);
    if (base!=MAP_FAILED && fileLen>0 && mmap(base, fileLen, PROT_READ, MAP_PRIVATE|MAP_FIXED, fd, 0)==MAP_FAILED) {
        munmap(base, end);
        base=(char*) MAP_FAILED;
    };
    close(fd);
    return base==MAP_FAILED ? NULL : base+offset;
};
