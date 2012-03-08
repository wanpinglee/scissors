
#ifndef SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_
#define SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_

#include <string>

using std::string;

struct Parameters {
  // i/o parameters
  string input_bam;             // -i  --input
  string input_reference_hash;  // -r  --reference-hash-table
  string output_bam;            // -o  --output

  string reference_filename;
  string hash_filename;

  // operation parameters
  int   fragment_length;  // -l --fragmenr-length
  int   mate_window_size; // -w --window-size
  float allowed_clip;     // -c --allowed-clip
  bool  is_input_sorted;  // -s --is-input-sorted
  int   processors;       // -p --processors
  bool  detect_special;   // -S --special-reference

  // alignment filter
  float aligned_base_rate;
  float allowed_mismatch_rate;
  int   trimming_match_score;
  int   trimming_mismatch_penalty;
  int   trimming_gap_penalty;
	
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
      , mate_window_size(2)
      , allowed_clip(0.2)
      , is_input_sorted(false)
      , processors(1)
      , detect_special(false)
      , aligned_base_rate(0.3)
      , allowed_mismatch_rate(0.1)
      , trimming_match_score(1)
      , trimming_mismatch_penalty(2)
      , trimming_gap_penalty(2)
  {}
};

void ParseArgumentsOrDie(const int argc, char* const * argv, Parameters* param);

#endif
