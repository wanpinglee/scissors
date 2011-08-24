#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <string>

#include "string_converter.h"
#include "parameter_parser.h"

using std::cout;
using std::endl;
using std::string;

ParameterParser::ParameterParser(const int argc, char* const * argv) {
	// parse the arguments and save the parameters
	ParseArgumentsOrDie( argc, argv );
	
	// check that all required parameters are set
	if ( !CheckParameters() )
		exit(1);
}

void ParameterParser::ParseArgumentsOrDie(const int argc, char* const * argv) {

	if ( argc == 1 ) { // no argument
		PrintHelp( argv );
		exit( 1 );
	}

	// record command line
	command_line = argv[0];
	for ( int i = 1; i < argc; ++i ) {
		command_line += " ";
		command_line += argv[i];
	}

	const char *short_option = "hi:r:o:l:s";

	const struct option long_option[] = {
		{ "help", no_argument, NULL, 'h' },
		{ "input", required_argument, NULL, 'i' },
		{ "output", required_argument, NULL, 'o' },
		{ "reference-hash-table", required_argument, NULL, 'r' },

		{ "fragment-length", required_argument, NULL, 'l' },
		{ "is-input-sorted", no_argument, NULL, 's'},

		{ 0, 0, 0, 0 }
	};

	int c = 0;
	bool help = false;
	while ( true ) {
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
				input_bam =  optarg;
				break;
			case 'o':
				output_bam = optarg;
				break;
			case 'r':
				input_reference_hash = optarg;
				reference_filename   = input_reference_hash + ".ref"; 
				hash_filename        = input_reference_hash + ".ht";
				break;
			// operation parameters
			case 'l':
				if ( !convert_from_string( optarg, fragment_length ) )
					cout << "WARNING: Cannot parse --fragment-length." << endl;
				break;
			case 's':
				is_input_sorted = true;
			default:
				break;
		}

	}

	if ( help ) {
		PrintHelp( argv );
		exit( 1 );
	}
}

// true: passing the checker
bool ParameterParser::CheckParameters(void) {
	
	bool errorFound = false;
	// necessary parameters
	if ( input_bam.empty() ) {
		cout << "ERROR: Please specific an input file, -i." << endl;
		errorFound = true;
	}
	
	if ( output_bam.empty() ) {
		cout << "ERROR: Please specific an output file, -o." << endl;
		errorFound = true;
	}
	
	if ( input_reference_hash.empty() ) {
		cout << "ERROR: Please specific a reference hash table, -r." << endl;
		errorFound = true;
	}

	return !errorFound;

}

void ParameterParser::PrintHelp(const char* const * argv) {
	cout
		<< endl
		<< "usage: " << argv[0] << " [OPTIONS] -i <FILE> -o <FILE> -r <FILE>"
		<< endl
		<< endl
		<< "Help:" << endl
		<< endl
		<< "   -h --help                 Print this help dialog." << endl
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
		<< "                         Fragment length."
		<< "   -s --is-input-sorted" << endl

		<< endl;
}

