#ifndef DATASTRUCTURES_BAMALIGNMENT_H_
#define DATASTRUCTURES_BAMALIGNMENT_H_

#include <string>

#include <stdint.h>

using namespace std;

const unsigned short kBamCoreSize   = 32;
// cigar operator for bam packed cigar
const unsigned short kBamCmatch     = 0;
const unsigned short kBamCins       = 1;
const unsigned short kBamCdel       = 2;
const unsigned short kBamCrefSkip   = 3;
const unsigned short kBamCsoftClip  = 4;
const unsigned short kBamChardClip  = 5;
const unsigned short kBamCpad       = 6;

struct BamAlignment {
	string        query_name;
	unsigned char flag;
	int32_t       reference_index;
	int32_t       reference_begin;  // position
	int32_t       reference_end;    // used for bin calculation
	unsigned char mapping_quality;
	string        bam_packed_cigar;
	int32_t       mate_reference_index;
	int32_t       mate_position;
	int32_t       read_length;
	int32_t       isize;            // can be obtained by position & matePosition
	string        encoded_sequence;
	string        qual;
	string        tag;

	BamAlignment( void )
		: flag(0)
		, reference_index(-1)
		, reference_begin(-1)
		, reference_end(0)
		, mapping_quality(0)
		, mate_reference_index(-1)
		, mate_position(-1)
		, read_length(0)
		, isize(0)
	{}
};

#endif
