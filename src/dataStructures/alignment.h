#ifndef DATASTRUCTURES_ALIGNMENT_H_
#define DATASTRUCTURES_ALIGNMENT_H_

#include <string>
#include <stdint.h>

using std::string;

#define MOSAIK_NUM_NUCLEOTIDES 26

namespace Scissors {
namespace AlignmentConstant {

const char GAP = '-';
const float FLOAT_NEGATIVE_INFINITY = (float)-1e+30;

} // namespace Alignment


struct Alignment{
	string reference;
	string query;
	uint32_t reference_begin;
	uint32_t reference_end;
	uint32_t query_begin;
	uint32_t query_end;
	uint32_t num_mismatches;
	unsigned char quality;

	bool is_seq_inverse;
	bool is_seq_complement;

	Alignment()
		: reference_begin(0)
		, reference_end(0)
		, query_begin(0)
		, query_end(0)
		, num_mismatches(0)
		, quality(0)

		, is_seq_inverse(false)
		, is_seq_complement(false)
	{}
}; // Alignment
} //namespace Scissors
#endif // DATASTRUCTURES_ALIGNMENT_H_
