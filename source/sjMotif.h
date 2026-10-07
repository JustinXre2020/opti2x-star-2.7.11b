#ifndef H_sjMotif
#define H_sjMotif

#include "IncludeDefine.h"

// Splice junction motif from the first two (donor d1 d2) and last two (acceptor a1 a2) intron bases, coded as in
// Transcript::canonSJ: 1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT, 0 non-canonical (incl. any non-ACGT base).
// A table lookup replaces the motif if-chain in the stitching junction scan.
class SjMotifTable {
public:
    SjMotifTable();
    static uint8 base(char b) {//collapse all non-ACGT genome codes into 4
        return (uint8)b<4 ? (uint8)b : 4;
    };
    uint8 operator()(uint8 d1, uint8 d2, uint8 a1, uint8 a2) const {//arguments from base()
        return code[((d1*5+d2)*5+a1)*5+a2];
    };
private:
    uint8 code[5*5*5*5];
};

extern const SjMotifTable sjMotifTable;

#endif
