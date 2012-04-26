#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <string>

#include "string_converter.h"
#include "parameter_parser.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

bool CheckParameters(Parameters* param);
void PrintLongHelp(const string& program);
void PrintBriefHelp(const string& program);

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

	const char *short_option = "hi:r:o:l:w:c:sp:SQ:B:M:";

	const struct option long_option[] = {
		// long help
		{"help", no_argument, NULL, 1},

		// i/o parameters
		{"input", required_argument, NULL, 'i'},
		{"output", required_argument, NULL, 'o'},
		{"reference-hash-table", required_argument, NULL, 'r'},

		// operation parameters
		{"fragment-length", required_argument, NULL, 'l'},
		{"window-size", required_argument, NULL, 'w'},
		{"allowed-clip", required_argument, NULL, 'c'},
		{"is-input-sorted", no_argument, NULL, 's'},
		{"processors", required_argument, NULL, 'p'},
		{"special-insertion", no_argument, NULL, 'S'},
		{"not-medium-sized-indel", no_argument, NULL, 5},

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
			case 'r':
				param->input_reference_hash = optarg;
				param->reference_filename   = 
				  param->input_reference_hash + ".ref"; 
				param->hash_filename        = 
				  param->input_reference_hash + ".ht";
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
			case 's':
				param->is_input_sorted = true;
			case 'p':
				if (!convert_from_string(optarg, param->processors))
					cerr << "WARNING: Cannot parse -p --processors." << endl;
				break;
			case 'S':
				param->detect_special = true;
				break;
			case 5:
				param->not_medium_sized_indel = true;

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
		cerr << "ERROR: Please specific an input file, -i." << endl;
		errorFound = true;
	}
	
	if (param->output_bam.empty()) {
		cerr << "ERROR: Please specific an output file, -o." << endl;
		errorFound = true;
	}
	
	if (param->input_reference_hash.empty()) {
		cerr << "ERROR: Please specific a reference hash table, -r." << endl;
		errorFound = true;
	}

	if (param->fragment_length < 1) {
		cerr << "ERROR: Please specific the fragment length, -l." << endl
		     << "       The value should be greater than 0." << endl;
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

void PrintBriefHelp(const string& program) {
	cout
		<< endl
		<< "usage: " << program << " [OPTIONS] -i <FILE> -o <FILE> -r <FILE_PREFIX> -l <INT>"
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
		<< "   -r --reference-hash-table <FILE_PREFIX>" << endl 
		<< "                         Hash table of the genome. A reference file (FILE_PREFI-" << endl
		<< "                         X.ref) and a hash table (FILE_PREFIX.ht) will be loaded" << endl
		<< "                         ." << endl
		<< endl

		<< "Operations:" << endl
		<< endl
		<< "   -l --fragment-length <INT>" << endl
		<< "                         Fragment length." << endl
		<< endl;
		

}

void PrintLongHelp(const string& program) {
	cout
		<< endl
		<< "usage: " << program << " [OPTIONS] -i <FILE> -o <FILE> -r <FILE_PREFIX> -l <INT>"
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
		<< "   -r --reference-hash-table <FILE_PREFIX>" << endl 
		<< "                         Hash table of the genome. A reference file (FILE_PREFI-" << endl
		<< "                         X.ref) and a hash table (FILE_PREFIX.ht) will be loaded" << endl
		<< "                         ." << endl
		<< endl
		
		<< "Operations:" << endl
		<< endl
		<< "   -l --fragment-length <INT>" << endl
		<< "                         Fragment length." << endl
		<< "   -w --window-size <INT>" << endl
		<< "                         Window size for searching mates. [fragment_length * 2]" << endl
		<< "   --discovery-window-size <INT>" << endl
		<< "                         Window size for discovering events. [10000]" << endl
		<< "   -s --is-input-sorted" << endl
		<< "   -p --processors <INT> Use # of processors." << endl
		<< "   -S --special-insertion" << endl
		<< "                         Detect insertions in special references, e.g. MEIs." << endl
		<< endl

		<< "Original BAM alignments filters:" << endl
		<< endl
		<< "   -Q --mapping-quality-threshold <INT>" << endl
		<< "                         Mapping quality threshold (0 - 255) of candidatei alig-" << endl
		<< "                         nments. [10]" << endl
		<< "   -c --allowed-clip <FLOAT>"  << endl
		<< "                         Percentage (0.0 - 1.0) of allowed soft clip of anchors." << endl
		<< "                         [0.2]" << endl
		<< endl

		<< "Split-read alignment filters:" << endl
		<< endl
		<< "   -B --aligned-base-rate <FLOAT>" << endl
		<< "                         Minimum aligned-base rate (0.0 - 1.0) of split-read al-" << endl
		<< "                         ignments. [0.3]" << endl
		<< "   -M --allowed-mismatch-rate <FLOAT>" << endl
		<< "                         Maximum mismatch rate (0.0 - 1.0) in split-read alignm-" << endl
		<< "                         ents. [0.1]" << endl
		<< "   --trimming-match-score <INT>" << endl
		<< "                         Match score for alignment trimming. [1]" << endl
		<< "   --trimming-match-penalty <INT>" << endl
		<< "                         Mismatch penalty for alignment trimming. [2]" << endl
		<< "   --trimming-gap-penalty <INT>" << endl
		<< "                         Gap penalty for alignment trimming. [2]" << endl

		<< endl;
}

