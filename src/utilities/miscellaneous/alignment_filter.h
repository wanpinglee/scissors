#ifndef UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
#define UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_

struct Alignment;

struct AlignmentFilter {
  int trimming_match_score;
  int trimming_mismatch_score;
  int trimming_gap_score;
};

namespace AlignmentFilterApplication {
bool TrimAlignment(const AlignmentFilter& filter, Alignment* al);
} // AlignmentFilterApplication

#endif // UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
