#ifndef UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
#define UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_

struct Alignment;

namespace AlignmentFilter {
struct Filter {
  int trimming_match_score;
  int trimming_mismatch_score;
  int trimming_gap_score;
};

bool TrimAlignment(const int& match_score, const int& mismatch_score, 
                   const int& gap_score, Alignment* al);

} // namespace AlignmentFilter
#endif // UTILITY_MISCELLANEOUS_ALIGNMENT_FILTER_H_
