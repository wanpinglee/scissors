
#ifndef _ParameterParser_H_
#define _ParameterParser_H_

#include <string>

using namespace std;

class ParameterParser {
	
	public:
	// i/o parameters
	string inputBam;            // -i  --input
	string inputReferenceHash;  // -r  --reference-hash-table
	string outputBam;           // -o  --output

	// operation parameters
	
	// command line
	const string commandLine;
	

	// functions
	ParameterParser ( const int argc, char* const * argv );
	void printHelp ( const char* const * argv );
	bool checkParameters ( void );
};

#endif
