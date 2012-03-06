#include "dataStructures/alignment.h"

namespace AlignmentFilter {
bool TrimAlignment(
    const int& match_score,
    const int& mismatch_score,
    const int& gap_score,
    Alignment* al){

  // find the sublist whose scores is max
  int maxi = 0, max = 0, start = 0, length = 0, best_start = 0, best_length = 0;
  for (unsigned int i = 0; i < al->reference.size(); ++i) {
    // get score
    int score = 0;
    if ((al->reference[i] == AlignmentConstant::GAP) || (al->query[i] == AlignmentConstant::GAP)) {
      score = gap_score;
    } else {
      if (al->reference[i] != al->query[i]) score = mismatch_score;
      else score = match_score;
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
  if (length > 0) {
    int end = start + length - 1;
    //sanity check
    if (end > (al->reference.size() - 1)) {
      return false;
    } else {
      if (end < (al->reference.size() - 1)) {
        for (unsigned int i = end + 1; i < al->reference.size(); ++i) {
	  if (al->reference[i] != AlignmentConstant::GAP)
	    --(al->reference_end);
	  if (al->query[i] != AlignmentConstant::GAP)
	    --(al->query_end);
	}
	int size = al->reference.size();
	al->reference.erase(end + 1, size - end - 1);
	al->query.erase(end + 1, size - end - 1);
      }
      
      if (start > 0) {
        for (int i = 0; i < start; ++i) {
	  if (al->reference[i] != AlignmentConstant::GAP)
	    ++(al->reference_begin);
	  if (al->query[i] != AlignmentConstant::GAP)
	    ++(al->query_begin);
	}
        al->reference.erase(0, start);
        al->query.erase(0, start);
      }
      return true;
    }
  } else {
    return false;
  }
} // TrimAlignment

} // namespace AlignmentFilter
