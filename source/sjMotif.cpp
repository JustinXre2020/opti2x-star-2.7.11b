#include "sjMotif.h"
#include <cstring>

const SjMotifTable sjMotifTable;

SjMotifTable::SjMotifTable() {
    memset(code, 0, sizeof(code));
    const uint8 motifs[6][4]={{2,3,0,2},  //1 GTAG
                              {1,3,0,1},  //2 CTAC
                              {2,1,0,2},  //3 GCAG
                              {1,3,2,1},  //4 CTGC
                              {0,3,0,1},  //5 ATAC
                              {2,3,0,3}}; //6 GTAT
    for (uint8 im=0; im<6; im++) {
        const uint8 *m=motifs[im];
        code[((m[0]*5+m[1])*5+m[2])*5+m[3]]=im+1;
    };
};
