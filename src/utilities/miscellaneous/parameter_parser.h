#ifndef SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_
#define SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_

#include <string>
#include "dataStructures/technology.h"

using std::string;

namespace Scissors {
struct Parameters {
  // i/o parameters
  string input_bam;             // -i  --input
  string input_reference_fasta; // -f  --fasta
  string input_special_fasta;   // -s  --special-fasta
  string output_bam;            // -o  --output
  string output_complete_bam;   // -O  --complete-bam

  // operation parameters
  int   fragment_length;        // -l --fragmenr-length
  int   mate_window_size;       // -w --window-size; default: fragment_length * 2
  int   discovery_window_size;  // --discovery-window-size
  bool  is_input_sorted;        // --is-input-sorted
                                // getopt returns 6
  int   processors;             // -p --processors
  bool  detect_special;         // when -s <FASTA> is given 
  bool  not_medium_sized_indel; // --not-medium-sized-indel
                                // getopt returns 5
  bool  not_special_insertion_inversion; // --not-special-insertion-inversion
                                         // getopt returns 7
  Technology technology;        // -t --technology

  // original alignment filters
  int mapping_quality_threshold; // -Q --mapping-quality-threshold
  float allowed_clip;            // -c --allowed-clip
  string region;                 // -r --region

  // split-read alignment filters
  float aligned_base_rate;         // -B --aligned-base-rate
  float allowed_mismatch_rate;     // -M --allowed-mismatch-rate
  int   trimming_match_score;      // --trimming-match-score
                                   // getopt returns 2
  int   trimming_mismatch_penalty; // --trimming-mismatch-penalty
                                   // getopt returns 3
  int   trimming_gap_penalty;      // --trimming-gap-penalty
                                   // getopt returns 4
	
  // command line
  string command_line;

  // default values
  Parameters()
      : input_bam()
      , input_reference_fasta()
      , input_special_fasta()
      , output_bam()
      , output_complete_bam()
      , fragment_length(-1)
      , mate_window_size(-1)
      , discovery_window_size(10000)
      , is_input_sorted(false)
      , processors(1)
      , detect_special(false)
      , not_medium_sized_indel(false)
      , not_special_insertion_inversion(false)
      , technology(TECH_NONE)
      , mapping_quality_threshold(10)
      , allowed_clip(0.2)
      , region()
      , aligned_base_rate(0.2)
      , allowed_mismatch_rate(0.1)
      , trimming_match_score(1)
      , trimming_mismatch_penalty(7)
      , trimming_gap_penalty(7)
      , command_line()
  {}
};

void ParseArgumentsOrDie(const int argc, char* const * argv, Parameters* param);
} // namespace Scissors
#endif
