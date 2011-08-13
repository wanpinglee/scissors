#ifndef _BAMHEADER_H_
#define _BAMHEADER_H_ 


const unsigned short BAM_CORE_SIZE  = 32;
// cigar operator for bam packed cigar
const unsigned short BAM_CMATCH     = 0;
const unsigned short BAM_CINS       = 1;
const unsigned short BAM_CDEL       = 2;
const unsigned short BAM_CREF_SKIP  = 3;
const unsigned short BAM_CSOFT_CLIP = 4;
const unsigned short BAM_CHARD_CLIP = 5;
const unsigned short BAM_CPAD       = 6;


#endif
