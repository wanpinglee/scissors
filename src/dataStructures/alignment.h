
#ifndef DATASTRUCTURES_ALIGNMENT_H_
#define DATASTRUCTURES_ALIGNMENT_H_

#include <string>

#include <stdint.h>

using std::string;

#define MOSAIK_NUM_NUCLEOTIDES 26

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

	bool TrimAlignment(
	    const int& match_score = 1,
	    const int& mismatch_score = -1,
	    const int& gap_score = -2)
	{
	  // find the sublist whose scores is max
	  int maxi = 0, max = 0, start = 0, length = 0;
	  for (unsigned int i = 0; i < reference.size(); ++i) {
	    // get score
	    int score = 0;
	    if ((reference[i] == AlignmentConstant::GAP) || (query[i] == AlignmentConstant::GAP)) {
	      score = match_score;
	    } else {
	      if (reference[i] != query[i]) score = mismatch_score;
	      else score = match_score;
	    }

	    // find max
	    if ((maxi + score) > 0) {
	      maxi += score;
	    } else {
	      maxi = 0;
	      start = i;
	      length = 0;
	    }
	    if (maxi > max) {
	      max = maxi;
	      length = i - start + 1;
	    }
	  }

	  // trim alignment
	  if (length > 0) {
	    int end = start + length - 1;
	    //sanity check
	    if (end > (reference.size() - 1)) {
	      return false;
	    } else {
	      int size = reference.size();
	      reference.erase(end + 1, size - end - 1);
	      query.erase(end + 1, size - end - 1);
	      reference_end -= (size - end - 1);
	      query_end -= (size - end - 1);
	      
	      reference.erase(0, start);
	      query.erase(0, start);
	      reference_begin += start;
	      query_begin += start;
	      return true;
	    }
	  }
	} // TrimAlignment
}; // Alignment

#endif // DATASTRUCTURES_ALIGNMENT_H_
