#ifndef H_fileMapping
#define H_fileMapping

#include "IncludeDefine.h"

// Maps length bytes of a file from offset read-only, private to this process; bytes past the end of the file read as
// zero. Returns NULL if the file cannot be opened or mapped. The mapping lives until the process exits.
char* mapFileReadOnly(const string& path, uint64 offset, uint64 length);

#endif
