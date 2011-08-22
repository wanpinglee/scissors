
#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
//#include "hasher/reader/SR_InHashTable.h"
//#include "hasher/reader/SR_Reference.h"
//#include "hasher/reader/SR_HashRegionTable.h"
#include "utilities/SR_BamInStream.h"
#include "hasher/common/SR_Types.h"
#include "dataStructures/SR_QueryRegion.h"
}

#include "utilities/bam_writer.h"
#include "utilities/parameter_parser.h"

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

	// Parse the arguments and store them
	// The program will exit(1) if any errors or
	//     missing required parameters are found
	const ParameterParser parameter_parser( argc, argv );
	
	// bam file reader
	SR_BamInStream* bam_reader = SR_BamInStreamAlloc( parameter_parser.input_bam.c_str() );

	// bam file writer
	BamWriter bam_writer( parameter_parser.output_bam );
	bam_writer.Open();

	// bam records are in SR_QueryRegion structure
	SR_QueryRegion* query_region = SR_QueryRegionAlloc();

	while( SR_BamInStreamGetPair( &(query_region->pAnchor), &(query_region->pOrphan), bam_reader ) == SR_OK ) {
		;
	}


	// free memory and close files
	SR_QueryRegionFree( query_region );
	SR_BamInStreamFree( bam_reader );
	bam_writer.Close();

	cout << "Program done." << endl;

	return 0;

}
