
#ifndef SRC_UTILITIES_BAMUTILITIES_H_
#define SRC_UTILITIES_BAMUTILITIES_H_

#include <stdint.h>

#include <string>

using std::string;

namespace BamUtilities {

// Given sequence, generate bam-format encoded sequence 
// and store it in encodedSequence
void EncodeQuerySequence( string& encodedSequence, const string& sequence );


// Given aligned reference and query sequences,
// generate bam-foramt packed cigar and store it in packed_cigar
bool GetPackedCigar( 
	string& packed_cigar, 
	const string& reference,
	const string& query,
	const uint32_t& query_begin,
	const uint32_t& query_end,
	const uint32_t& read_length);

} // namespace BamUtilities

#endif // SRC_UTILITIES_BAMUTILITIES_H_
