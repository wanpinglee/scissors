
#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
#include "hasher/reader/SR_InHashTable.h"
#include "hasher/reader/SR_Reference.h"
}

#include "utilities/ParameterParser.h"

using namespace std;

inline FILE* openFile( const string filename, const string type ) {

	FILE* file = fopen( filename.c_str(), type.c_str() );
	if ( file == NULL ) {
		cout << "ERROR: Cannot open the file " << filename << endl;
		exit( 1 );
	}
	
	return file;
}

int main ( int argc, char** argv ) {

	// parse the parameters and store them
	// the program will exit(1) if any errors are found
	const ParameterParser parameterparser( argc, argv );

	uint64_t referenceBufferLength = 1000000;
	SR_Reference* pRef = SR_ReferenceAlloc( referenceBufferLength );

	// open hash table
	FILE* hashTableInput = openFile( parameterparser.inputReferenceHash + ".ref", "rb" );

	unsigned short hashSize = SR_ReadHashSize( hashTableInput );

	// create a hash region table object
	// the parameter given here is the edge tolerance percentage which is not useful here
	HashRegionTable* pRegionTable = HashRegionTableAlloc(0.25);
	
	
	// close files
	fclose( hashTableInput );

	return 0;

}
