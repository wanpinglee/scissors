
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
	
  // command line
  string command_line;

  // default values
  Parameters()
      : mate_window_size(2)
      , allowed_clip(0.2)
      , is_input_sorted(false)
      , processors(1)
  {}
};

void ParseArgumentsOrDie(const int argc, char* const * argv, Parameters* param);

#endif