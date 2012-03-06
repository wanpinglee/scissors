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
void PrintHelp(const string& program);

void ParseArgumentsOrDie(const int argc, char* const * argv, 
    Parameters* param) {

	if (argc == 1) { // no argument
		PrintHelp(argv[0]);
		exit(1);
	}

	// record command line
	param->command_line = argv[0];
	for ( int i = 1; i < argc; ++i ) {
		param->command_line += " ";
		param->command_line += argv[i];
	}

	const char *short_option = "hi:r:o:l:w:c:sp:SCM";

	const struct option long_option[] = {
		{ "help", no_argument, NULL, 'h' },
		{ "input", required_argument, NULL, 'i' },
		{ "output", required_argument, NULL, 'o' },
		{ "reference-hash-table", required_argument, NULL, 'r' },

		{ "fragment-length", required_argument, NULL, 'l' },
		{ "window-size", required_argument, NULL, 'w' },
		{ "allowed-clip", required_argument, NULL, 'c' },
		{ "is-input-sorted", no_argument, NULL, 's'},
		{ "processors", required_argument, NULL, 'p'},
		{ "special", no_argument, NULL, 'S'},

		{ "alignment-coverage-rate", no_argument, NULL, 'C'},
		{ "allowed-mismatch-rate", no_argument, NULL, 'M'},

		{ 0, 0, 0, 0 }
	};

	int c = 0;
	bool help = false;
	while (true) {
		int optionIndex = 0;
		c = getopt_long( argc, argv, short_option, long_option, &optionIndex );
		
		if ( c == -1 ) // end of options
			break;
		
		switch ( c ) {
			// help
			case 'h':
				help = true;
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

			case 'C':
				if (!convert_from_string(optarg, param->alignment_coverage_rate))
					cerr << "WARNING: Cannot parse -C --alignment-coverage-rate." << endl;
				break;
			case 'M':
				if (!convert_from_string(optarg, param->allowed_mismatch_rate))
					cerr << "WARNING: Cannot parse -M --allowed-mismatch-rate." << endl;
				break;
			default:
				break;
		}

	}

	if (help) {
		PrintHelp(argv[0]);
		exit(1);
	}

	CheckParameters(param);
}

// true: passing the checker
bool CheckParameters(Parameters* param) {
	
	bool errorFound = false;
	// necessary parameters
	if ( param->input_bam.empty() ) {
		cerr << "ERROR: Please specific an input file, -i." << endl;
		errorFound = true;
	}
	
	if ( param->output_bam.empty() ) {
		cerr << "ERROR: Please specific an output file, -o." << endl;
		errorFound = true;
	}
	
	if ( param->input_reference_hash.empty() ) {
		cerr << "ERROR: Please specific a reference hash table, -r." << endl;
		errorFound = true;
	}

	if ( param->fragment_length == -1 ) {
		cerr << "ERROR: Please specific fragment length, -l." << endl;
		errorFound = true;
	}

	// unnecessary parameters
	if ( ( param->allowed_clip < 0.0 ) || ( param->allowed_clip > 1.0 ) ) {
		cerr << "WARNING: -c should be in [0.0 - 1.0]. Set it to default, 0.2." << endl;
		param->allowed_clip = 0.2;
	}

	if ( param->mate_window_size == 0 ) {
		cerr << "WARNING: -w should not be zero. Set it to default, 2." << endl;
		param->mate_window_size = 2;
	}

	return !errorFound;

}

void PrintHelp(const string& program) {
	cout
		<< endl
		<< "usage: " << program << " [OPTIONS] -i <FILE> -o <FILE> -r <FILE_PREFIX> -l <INT>"
		<< endl
		<< endl
		<< "Help:" << endl
		<< endl
		<< "   -h --help             Print this help dialog." << endl
		<< endl
		<< "Input & Output:" << endl
		<< endl
		<< "   -i --input <FILE>     Input BAM-format file." << endl
		<< "   -o --output <FILE>    Output BAM-format file." << endl
		<< "   -r --reference-hash-table <FILE>" << endl 
		<< "                         Hash table of the genome. A reference file (FILE.ref)" << endl
		<< "                         and a hash table (FILE.ht) will be loaded." << endl
		<< endl
		
		<< "Operations" << endl
		<< endl
		<< "   -l --fragment-length <INT>" << endl
		<< "                         Fragment length." << endl
		<< "   -w --window-size <INT>" << endl
		<< "                         Window size (-w x -l) for searching mates in bam." << endl
		<< "                         Default: 2" << endl
		<< "   -c --allowed-clip <FLOAT>"  << endl
		<< "                         Percentage [0.0 - 1.0] of allowed soft clip of anchors." << endl
		<< "                         Default: 0.2" << endl
		<< "   -s --is-input-sorted" << endl
		<< "   -p --processors <INT> Use # of processors." << endl
		<< "   -S --special-reference" << endl
		<< "                         Detect special references, e.g. MEI." << endl

		<< "Alignment filter" << endl
		<< endl
		<< "   -C --alignment-coverage-rate <FLOAT>" << endl
		<< "                         Minimum coverage rate of original reads covered by " << endl
		<< "                         split-read alignments." << endl
		<< "                         Default: 0.9" << endl
		<< "   -M --allowed-mismatch-rate <FLOAT>" << endl
		<< "                         Maximum mismatch rate in split-read alignments." << endl
		<< "                         Default: 0.1"

		<< endl;
}

