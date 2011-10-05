
#ifndef SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_
#define SRC_UTILITIES_MISCELLANEOUS_ParameterParser_H_

#include <string>

using std::string;

class ParameterParser {
 public:
  ParameterParser(const int argc, char* const * argv);

  // i/o parameters
  string input_bam;             // -i  --input
  string input_reference_hash;  // -r  --reference-hash-table
  string output_bam;            // -o  --output

  string reference_filename;
  string hash_filename;

  // operation parameters
  unsigned int fragment_length;  // -l --fragmenr-length
  unsigned int mate_window_size; // -w --window-size
  float allowed_clip;            // -c --allowed-clip
  bool is_input_sorted;          // -s --is-input-sorted
	
  // command line
  string command_line;
	

  // functions
  void PrintHelp(const char* const * argv);
	
 private:
  void ParseArgumentsOrDie(const int argc, char* const * argv);
  bool CheckParameters(void);
};

#endif
