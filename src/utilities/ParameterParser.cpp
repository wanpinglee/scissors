#include <stdlib.h>

#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "ParameterParser.h"

ParameterParser::ParameterParser
	( const int argc
	, char* const * argv ) {

	if ( argc == 1 ) { // no argument
		printHelp( argv );
		exit( 1 );
	}

	for ( int i = 0; i <= argc; i++ )
		commandLine += *argv[i];

	const char *short_option = "hi:r:o:";

	const struct option long_option[] = {
		{ "help", no_argument, NULL, 'h' },
		{ "input", required_argument, NULL, 'i'},
		{ "output", required_argument, NULL, 'o' },
		{ "reference-hash-table", required_argument, NULL, 'r'},

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
				inputBam =  optarg;
				break;
			case 'o':
				outputBam = optarg;
				break;
			case 'r':
				inputReferenceHash = optarg;
				break;
			default:
				break;
		}

	}

	if ( help ) {
		printHelp( argv );
		exit( 1 );
	}
	else if ( !checkParameters() )
		exit( 1 );
}

// true: passing the checker
bool ParameterParser::checkParameters( void ) {
	
	bool errorFound = false;
	// necessary parameters
	if ( inputBam.empty() ) {
		cout << "ERROR: Please specific an input file." << endl;
		errorFound = true;
	}
	
	if ( outputBam.empty() ) {
		cout << "ERROR: Please specific an output file." << endl;
		errorFound = true;
	}
	
	if ( inputReferenceHash.empty() ) {
		cout << "ERROR: Please specific a reference hash table." << endl;
		errorFound = true;
	}

	return errorFound;

}

void ParameterParser::printHelp( const char* const * argv ) {
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

		<< endl;
}

