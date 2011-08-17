
#ifndef DATASTRUCTURES_ALIGNMENT_H_
#define DATASTRUCTURES_ALIGNMENT_H_

#include <string>

#include <stdint.h>

using namespace std;

struct Alignment{
	string reference;
	string query;
	int32_t reference_begin;
	int32_t reference_end;
	int32_t query_begin;
	int32_t query_end;

	Alignment()
		: reference_begin(0)
		, reference_end(0)
		, query_begin(0)
		, query_end(0)
	{}
};

#endif // DATASTRUCTURES_ALIGNMENT_H_
