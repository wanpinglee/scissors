
#ifndef SRC_UTILITIES_BAMUTILITIES_H_
#define SRC_UTILITIES_BAMUTILITIES_H_

#include <stdint.h>

#include <string>
#include <vector>

#include "samtools/bam.h"

using std::string;
using std::vector;

namespace BamUtilities {

void ReplaceHeaderText( const bam_header_t* header );

// Given sequence, generate bam-format encoded sequence 
// and store it in encodedSequence
void EncodeQuerySequence( string& encodedSequence, const string& sequence );


// Given the packed_cigar,
// convert packed_cigar into uint8_t for writing
//inline const uint8_t* ConvertPackedCigar( const vector<uint32_t>& packed_cigar );

// Given aligned reference and query sequences,
// generate bam-foramt packed cigar and store it in packed_cigar
bool GetPackedCigar( 
	vector<uint32_t>& packed_cigar, 
	const string& reference,
	const string& query,
	const uint32_t& query_begin,
	const uint32_t& query_end,
	const uint32_t& read_length);

} // namespace BamUtilities

#endif // SRC_UTILITIES_BAMUTILITIES_H_
