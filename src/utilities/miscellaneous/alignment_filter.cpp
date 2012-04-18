#include "alignment_filter.h"

#include <iostream>
#include "dataStructures/alignment.h"

namespace AlignmentFilterApplication {
bool TrimAlignment(
    const AlignmentFilter& filter,
    Alignment* al){

  // find the sublist whose scores is max
  int maxi = 0, max = 0, start = 0, length = 0, best_start = 0, best_length = 0;
  for (unsigned int i = 0; i < al->reference.size(); ++i) {
    // get score
    int score = 0;
    if ((al->reference[i] == AlignmentConstant::GAP) || (al->query[i] == AlignmentConstant::GAP)) {
      score = filter.trimming_gap_score;
    } else {
      if (al->reference[i] != al->query[i]) score = filter.trimming_mismatch_score;
      else score = filter.trimming_match_score;
    }

    // find max
    if ((maxi + score) > 0) {
      maxi += score;
      ++length;
    } else {
      maxi = 0;
      start = i + 1;
      length = 0;
    }
    if (maxi > max) {
      max = maxi;
      best_start = start;
      best_length = length;
    }
  }

  // trim alignment
  if ((best_length > 0) 
  && (static_cast<unsigned int>(best_length) < al->reference.size())) {
    int end = best_start + best_length - 1;
    //sanity check
    if (end > (al->reference.size() - 1)) {
      return false;
    } else {
      // trim from the tail
      if (end < (al->reference.size() - 1)) {
        for (unsigned int i = end + 1; i < al->reference.size(); ++i) {
	  if (al->reference[i] != AlignmentConstant::GAP) --(al->reference_end);
	  if (al->query[i] != AlignmentConstant::GAP) --(al->query_end);
	  if (al->reference[i] != al->query[i]) --(al->num_mismatches);
	}
	int size = al->reference.size();
	al->reference.erase(end + 1, size - end - 1);
	al->query.erase(end + 1, size - end - 1);
      }
      
      // trim from the begin
      if (best_start > 0) {
        for (int i = 0; i < best_start; ++i) {
	  if (al->reference[i] != AlignmentConstant::GAP) ++(al->reference_begin);
	  if (al->query[i] != AlignmentConstant::GAP) ++(al->query_begin);
	  if (al->reference[i] != al->query[i]) --(al->num_mismatches);
	}
        al->reference.erase(0, best_start);
        al->query.erase(0, best_start);
      }
      return true;
    }
  } else {
    return false;
  }
} // TrimAlignment

} // namespace AlignmentFilter
