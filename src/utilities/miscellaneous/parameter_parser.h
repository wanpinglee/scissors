#ifndef SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_
#define SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_

#include <string>
#include "dataStructures/technology.h"

using std::string;

namespace Scissors {
struct Parameters {
  // i/o parameters
  string input_bam;             // -i  --input
  string input_reference_hash;  // -r  --reference-hash-table
  string output_bam;            // -o  --output

  string reference_filename;
  string hash_filename;

  // operation parameters
  int   fragment_length;        // -l --fragmenr-length
  int   mate_window_size;       // -w --window-size; default: fragment_length * 2
  int   discovery_window_size;  // --discovery-window-size
  bool  is_input_sorted;        // -s --is-input-sorted
  int   processors;             // -p --processors
  bool  detect_special;         // -S --special-reference
  bool  not_medium_sized_indel; // --not-medium-sized-indel
  Technology technology;        // -t --technology

  // original alignment filters
  int mapping_quality_threshold; // -Q --mapping-quality-threshold
  float allowed_clip;            // -c --allowed-clip
  string region;                 // -R --region

  // split-read alignment filters
  float aligned_base_rate;         // -B --aligned-base-rate
  float allowed_mismatch_rate;     // -M --allowed-mismatch-rate
  int   trimming_match_score;      // --trimming-match-score
  int   trimming_mismatch_penalty; // --trimming-mismatch-penalty
  int   trimming_gap_penalty;      // --trimming-gap-penalty
	
  // command line
  string command_line;

  // default values
  Parameters()
      : input_bam()
      , input_reference_hash()
      , output_bam()
      , reference_filename()
      , hash_filename()
      , fragment_length(-1)
      , mate_window_size(-1)
      , discovery_window_size(10000)
      , is_input_sorted(false)
      , processors(1)
      , detect_special(false)
      , not_medium_sized_indel(false)
      , technology(TECH_NONE)
      , mapping_quality_threshold(10)
      , allowed_clip(0.2)
      , aligned_base_rate(0.3)
      , allowed_mismatch_rate(0.1)
      , trimming_match_score(1)
      , trimming_mismatch_penalty(2)
      , trimming_gap_penalty(2)
  {}
};

void ParseArgumentsOrDie(const int argc, char* const * argv, Parameters* param);
} // namespace Scissors
#endif
