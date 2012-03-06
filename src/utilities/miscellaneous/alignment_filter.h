#ifndef UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
#define UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_

#include "dataStructures/alignment.h"

struct AlignmentFilter {
  float alignment_coverage_rate;
  float allowed_mismatch_rate;
  int trimming_match_score;
  int trimming_mismatch_score;
  int trimming_gap_score;
};

namespace AlignmentFilterApplication {
bool TrimAlignment(const AlignmentFilter& filter, Alignment* al);

inline bool FilterByMismatch(const AlignmentFilter& filter, const Alignment& al);

bool FilterByCoverage(const AlignmentFilter& filter, const int& original_length,
                      const Alignment& al);
} // AlignmentFilterApplication

bool AlignmentFilterApplication::FilterByMismatch(const AlignmentFilter& filter, const Alignment& al) {
  float mismatches = al.num_mismatches * filter.allowed_mismatch_rate;
  if (al.num_mismatches > mismatches) return true;
  else return false;
}

#endif // UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
