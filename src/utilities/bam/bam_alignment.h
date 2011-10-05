#ifndef DATASTRUCTURES_BAMALIGNMENT_H_
#define DATASTRUCTURES_BAMALIGNMENT_H_

#include <string>
#include <vector>

#include <stdint.h>

using std::string;
using std::vector;

namespace BamAlignmentConstant {

const unsigned char kBamCoreSize   = 32;
const unsigned char kBamCigarShift = 4;

// cigar operator for bam packed cigar
const unsigned char kBamCmatch     = 0;
const unsigned char kBamCins       = 1;
const unsigned char kBamCdel       = 2;
const unsigned char kBamCrefSkip   = 3;
const unsigned char kBamCsoftClip  = 4;
const unsigned char kBamChardClip  = 5;
const unsigned char kBamCpad       = 6;

} // namespace BamAlignmentConstant

struct BamAlignment {
	string           query_name;
	unsigned char    flag;
	int32_t          reference_index;
	int32_t          reference_begin;  // position
	int32_t          reference_end;    // used for bin calculation
	unsigned char    mapping_quality;
	vector<uint32_t> bam_packed_cigar;
	int32_t          mate_reference_index;
	int32_t          mate_position;
	int32_t          read_length;
	int32_t          isize;            // can be obtained by position & matePosition
	string           encoded_sequence;
	string           qual;
	string           tag;

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
