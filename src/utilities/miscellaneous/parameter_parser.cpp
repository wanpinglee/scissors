#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>
#include <string>

#include "string_converter.h"
#include "parameter_parser.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace Scissors {
bool CheckParameters(Parameters* param);
void PrintLongHelp(const string& program);
void PrintBriefHelp(const string& program);
void Convert_Technology(const string& optarg, Technology* technology);

void ParseArgumentsOrDie(const int argc, char* const * argv, 
    Parameters* param) {

	if (argc == 1) { // no argument
		PrintBriefHelp(argv[0]);
		exit(1);
	}

	// record command line
	param->command_line = argv[0];
	for ( int i = 1; i < argc; ++i ) {
		param->command_line += " ";
		param->command_line += argv[i];
	}

	const char *short_option = "hi:f:s:o:O:l:w:c:r:p:Q:B:M:Pt:";

	const struct option long_option[] = {
		// long help
		{"help", no_argument, NULL, 1},

		// i/o parameters
		{"input", required_argument, NULL, 'i'},
		{"output", required_argument, NULL, 'o'},
		{"complete-bam", required_argument, NULL, 'O'},
		{"fasta", required_argument, NULL, 'f'},
		{"special-fasta", required_argument, NULL, 's'},

		// operation parameters
		{"fragment-length", required_argument, NULL, 'l'},
		{"window-size", required_argument, NULL, 'w'},
		{"allowed-clip", required_argument, NULL, 'c'},
		{"region", required_argument, NULL, 'r'},
		{"is-input-sorted", no_argument, NULL, 6},
		{"processors", required_argument, NULL, 'p'},
		{"use-poor-mapped-mate", no_argument, NULL, 'P'},
		{"not-medium-sized-indel", no_argument, NULL, 5},
		{"not-special-insertion-inversion", no_argument, NULL, 7},
		{"technology", required_argument, NULL, 't'},

		// original bam alignment filters
		{"mapping-quality-threshold", no_argument, NULL, 'Q'},

		// split-read alignment filters
		{"aligned-base-rate", no_argument, NULL, 'B'},
		{"allowed-mismatch-rate", no_argument, NULL, 'M'},
		{"trimming-match-score", required_argument, NULL, 2},
		{"trimming-mismatch-penalty", required_argument, NULL, 3},
		{"trimming-gap-penalty", required_argument, NULL, 4},


		{0, 0, 0, 0}
	};

	int c = 0;
	bool brief_help = false;
	bool long_help  = false;
	while (true) {
		int option_index = 0;
		c = getopt_long(argc, argv, short_option, long_option, &option_index);
		
		if (c == -1) // end of options
			break;
		
		switch (c) {
			case 0:
				break;
			// help
			case 'h':
				brief_help = true;
				break;
			case 1:
				brief_help = false;
				long_help = true;
				break;
			// i/o parameters
			case 'i':
				param->input_bam =  optarg;
				break;
			case 'o':
				param->output_bam = optarg;
				break;
			case 'O':
				param->output_complete_bam = optarg;
				break;
			case 'f':
				param->input_reference_fasta = optarg;
				break;
			case 's':
				param->input_special_fasta = optarg;
				param->detect_special = true;
				break;
			// operation parameters
			case 'l':
				if (!convert_from_string(optarg, param->fragment_length))
					cerr << "WARNING: Cannot parse -l --fragment-length." << endl;
				break;
			case 'w':
				if (!convert_from_string(optarg, param->mate_window_size))
					cerr << "WARNING: Cannot parse -w --window-size." << endl;
				break;
			case 'c':
				if (!convert_from_string(optarg, param->allowed_clip))
					cerr << "WARNING: Cannot parse -c --allowed-clip." << endl;
				break;
			case 'r':
				param->region = optarg;
				break;
			case 6:
				param->is_input_sorted = true;
			case 'p':
				if (!convert_from_string(optarg, param->processors))
					cerr << "WARNING: Cannot parse -p --processors." << endl;
				break;
			case 'P': param->use_poor_mapped_mate = true;
			        break;
			case 5:
				param->not_medium_sized_indel = true;
				break;
			case 7: param->not_special_insertion_inversion = false;
				break;
			case 't': 
			        Convert_Technology(optarg, &(param->technology));
			        break;

			// original bam alignment filters
			case 'Q':
				if (!convert_from_string(optarg, param->mapping_quality_threshold))
					cerr << "WARNING: Cannot parse -Q --mapping-quality-threshold." << endl;
				break;

			// split-read alignment filters
			case 'B':
				if (!convert_from_string(optarg, param->aligned_base_rate))
					cerr << "WARNING: Cannot parse -B --aligned-base-rate." << endl;
				break;
			case 'M':
				if (!convert_from_string(optarg, param->allowed_mismatch_rate))
					cerr << "WARNING: Cannot parse -M --allowed-mismatch-rate." << endl;
				break;

			case 2:
				if (!convert_from_string(optarg, param->trimming_match_score))
					cerr << "WARNING: Cannot parse the argument of --trimming-match-score." << endl;
				break;
			case 3:
				if (!convert_from_string(optarg, param->trimming_mismatch_penalty))
					cerr << "WARNING: Cannot parse the argument of --trimming-mismatch-penalty." << endl;
				break;
			case 4:
				if (!convert_from_string(optarg, param->trimming_gap_penalty))
					cerr << "WARNING: Cannot parse the argument of --trimming-gap-penalty." << endl;
				break;
			

			default:
				cerr << "WARNING: Unkonw parameter: " << long_option[option_index].name << endl;
				break;
		}

	}

	if ((param->fragment_length != -1) && (param->mate_window_size == -1))
	  param->mate_window_size = param->fragment_length * 2;

	if (long_help) {
		PrintLongHelp(argv[0]);
		exit(1);
	}
	if (brief_help) {
		PrintBriefHelp(argv[0]);
		exit(1);
	}

	if (!CheckParameters(param)) exit(1);
}

// true: passing the checker
bool CheckParameters(Parameters* param) {
	
  bool errorFound = false;
  // necessary parameters
  if (param->input_bam.empty()) {
    cerr << "ERROR: Please specify an input BAM file, -i." << endl;
    errorFound = true;
  }
	
  if (param->output_bam.empty()) {
    cerr << "ERROR: Please specify an output BAM file, -o." << endl;
    errorFound = true;
  }
	
  if (param->input_reference_fasta.empty()) {
    cerr << "ERROR: Please specify a fasta file of references, -f." << endl;
    errorFound = true;
  }

  if (param->fragment_length < 1) {
    cerr << "ERROR: Please specify the fragment length, -l." << endl
         << "       The value should be greater than 0." << endl;
    errorFound = true;
  }

  if (param->technology == TECH_NONE) {
    cerr << "ERROR: Please specify the technology, -t." << endl
         << "       It should be ILLUMINA, 454, or SOLID." << endl;
    errorFound = true;
  }

  if(param->use_poor_mapped_mate && param->output_complete_bam.empty()) {
    cerr << "ERROR: Please specify the complete bam, -o," << endl
         << "       since -b is enabled." << endl;
    errorFound = true;
  }

  // unnecessary parameters
  if ((param->allowed_clip < 0.0) || (param->allowed_clip > 1.0)) {
    cerr << "WARNING: -c should be in [0.0 - 1.0]. Set it to default, 0.2." << endl;
    param->allowed_clip = 0.2;
  }

  if ((param->mate_window_size < 1) && (param->fragment_length > 0)) {
    cerr << "WARNING: -w should be greater 0. Set it to default, fragment_length * 2." << endl;
    param->mate_window_size = param->fragment_length * 2;
  }

  if (param->discovery_window_size < 0) {
    cerr << "WARNING: --discovery-window-size should be greater 0." << endl
         << "         Set it to default, 1000." << endl;
  }

  if ((param->aligned_base_rate < 0.0) || (param->aligned_base_rate > 1.0)) {
    cerr << "WARNING: -B should be in [0.0 - 1.0]. Set it to default, 0.3." << endl;
    param->aligned_base_rate = 0.3;
  }

  if ((param->allowed_mismatch_rate < 0.0) || (param->allowed_mismatch_rate > 1.0)) {
    cerr << "WARNING: -M should be in [0.0 - 1.0]. Set it to default, 0.1." << endl;
    param->allowed_mismatch_rate = 0.1;
  }

  return !errorFound;

}

void Convert_Technology(const string& optarg, Technology* technology) {
  string tech;
  tech.resize(optarg.size());
  std::transform(optarg.begin(), optarg.end(),tech.begin(), ::toupper);

  if (tech == "ILLUMINA")
    *technology = TECH_ILLUMINA;
  else if (tech == "454") 
    *technology = TECH_454;
  else if (tech == "SOLID") 
    *technology = TECH_SOLID;
  else 
    *technology = TECH_NONE;
  
}

void PrintBriefHelp(const string& program) {
	cout
		<< endl
		<< "usage: " << program << " [OPTIONS] -i <FILE> -o <FILE> -f <FILE> -l <INT> -t <STR>"
		<< endl
		<< "--help: for the complete help dialog."
		<< endl
		<< endl
		<< "Help:" << endl
		<< "   -h                    Print this help dialog." << endl
		<< "   --help                Print the complete help dialog." << endl
		<< endl
		<< "Input & Output:" << endl
		<< endl
		<< "   -i --input <FILE>     Input BAM file." << endl
		<< "   -o --output <FILE>    Output BAM file." << endl
		<< "   -f --fasta <FILE>     Input FASTA file."    << endl 
		<< endl

		<< "Operations:" << endl
		<< endl
		<< "   -l --fragment-length <INT>" << endl
		<< "                         Fragment length." << endl
		<< "   -t --technology <STR> ILLUMINA, 454, or SOLID." << endl
		<< endl;
		

}

void PrintLongHelp(const string& program) {
	cout
		<< endl
		<< "usage: " << program << " [OPTIONS] -i <FILE> -o <FILE> -r <FILE_PREFIX> -l <INT> -t <STR>"
		<< endl
		<< endl
		<< "Help:" << endl
		<< endl
		<< "   -h                    Print the brief help dialog." << endl
		<< "   --help                Print this help dialog." << endl
		<< endl
		<< "Input & Output:" << endl
		<< endl
		<< "   -i --input <FILE>     Input BAM file." << endl
		<< "   -o --output <FILE>    Output BAM file." << endl
		<< "   -O --complete-bam <FILE>" << endl
		<< "                         A generated bam contains original records" << endl
		<< "                         and alignments rescued by this split-read aligner." << endl
		<< "   -f --fasta            Input FASTA file." << endl
		<< "   -s --special-fasta <FILE>" << endl
		<< "                         A FASTA file of insertion sequences." << endl
		<< "                         Detect insertions in special references, e.g. MEIs." << endl
		<< endl
		
		<< "Operations:" << endl
		<< endl
		<< "   -l --fragment-length <INT>" << endl
		<< "                         Fragment length." << endl
		<< "   -w --window-size <INT>" << endl
		<< "                         Window size for searching mates. [fragment_length * 2]" << endl
		<< "   --discovery-window-size <INT>" << endl
		<< "                         Window size for discovering events. [10000]" << endl
		<< "   --is-input-sorted" << endl
		<< "   -p --processors <INT> Use # of processors." << endl
		<< "   -P --use-poor-mapped-mate" << endl
		<< "                         Use pairs with one mare good and the other mate that" << endl
		<< "                         are mapped but cannot pass -Q and -c filters." << endl
		<< "   --not-medium-sized-indel" << endl
		<< "   --not-special-insertion-inversion" << endl
		<< "                         When -s is given, the default is on." << endl
		<< "   -t --technology <STR> ILLUMINA, 454, or SOLID." << endl
		<< endl

		<< "Original BAM alignments filters:" << endl
		<< endl
		<< "   -Q --mapping-quality-threshold <INT>" << endl
		<< "                         Mapping quality threshold (0 - 255) of candidatei alig-" << endl
		<< "                         nments. [10]" << endl
		<< "   -c --allowed-clip <FLOAT>"  << endl
		<< "                         Percentage (0.0 - 1.0) of allowed soft clip of anchors." << endl
		<< "                         [0.2]" << endl
		<< "   -r --region <STR>     Targeted region; example: -r 1:500000-600000." << endl
		<< endl

		<< "Split-read alignment filters:" << endl
		<< endl
		<< "   -B --aligned-base-rate <FLOAT>" << endl
		<< "                         Minimum aligned-base rate (0.0 - 1.0) of split-read al-" << endl
		<< "                         ignments. [0.2]" << endl
		<< "   -M --allowed-mismatch-rate <FLOAT>" << endl
		<< "                         Maximum mismatch rate (0.0 - 1.0) in split-read alignm-" << endl
		<< "                         ents. [0.1]" << endl
		<< "   --trimming-match-score <INT>" << endl
		<< "                         Match score for alignment trimming. [1]" << endl
		<< "   --trimming-mismatch-penalty <INT>" << endl
		<< "                         Mismatch penalty for alignment trimming. [7]" << endl
		<< "   --trimming-gap-penalty <INT>" << endl
		<< "                         Gap penalty for alignment trimming. [7]" << endl

		<< endl;
}
} // namespace Scissors
