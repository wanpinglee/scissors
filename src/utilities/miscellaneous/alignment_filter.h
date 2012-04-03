#ifndef UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
#define UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_

#include <math.h>

#include "dataStructures/alignment.h"

struct AlignmentFilter {
  float aligned_base_rate;
  float allowed_mismatch_rate;
  int trimming_match_score;
  int trimming_mismatch_score;
  int trimming_gap_score;

  AlignmentFilter()
      : aligned_base_rate(0.1)
      , allowed_mismatch_rate(0.1)
      , trimming_match_score(1)
      , trimming_mismatch_score(-2)
      , trimming_gap_score(-2)
  {}
};

namespace AlignmentFilterApplication {
bool TrimAlignment(const AlignmentFilter& filter, Alignment* al);

inline bool FilterByMismatch(const AlignmentFilter& filter, const Alignment& al);

inline bool FilterByAlignedBaseThreshold(const AlignmentFilter& filter,
                                         const Alignment& al,
					 const int& read_length);
} // AlignmentFilterApplication

inline bool AlignmentFilterApplication::FilterByMismatch(
    const AlignmentFilter& filter, 
    const Alignment& al) {
  float allowed_mismatches = ceil(al.reference.size() * filter.allowed_mismatch_rate);
  if (al.num_mismatches < static_cast<unsigned int>(allowed_mismatches)) return true;
  else return false;
}

inline bool AlignmentFilterApplication::FilterByAlignedBaseThreshold(
    const AlignmentFilter& filter,
    const Alignment& al,
    const int& read_length) {
  int aligned_base_threshold = floor(read_length * filter.aligned_base_rate);
  int aligned_bases = al.query_end - al.query_begin + 1;
  if (aligned_bases > aligned_base_threshold) return true;
  else return false;
}

#endif // UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
