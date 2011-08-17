
#ifndef DATASTRUCTURES_ALIGNMENT_H_
#define DATASTRUCTURES_ALIGNMENT_H_

#include <string>

#include <stdint.h>

using namespace std;

struct Alignment{
	string reference;
	string query;
	uint32_t reference_begin;
	uint32_t reference_end;
	uint32_t query_begin;
	uint32_t query_end;
	uint32_t num_mismatches;

	bool is_reverse_complement;

	Alignment()
		: reference_begin(0)
		, reference_end(0)
		, query_begin(0)
		, query_end(0)
		, num_mismatches(0)

		, is_reverse_complement(false)
	{}
};

#endif // DATASTRUCTURES_ALIGNMENT_H_
