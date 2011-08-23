
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
#include "samtools/bam.h"
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

void Deconstruct( 
	SR_BamInStream* bam_reader,
	BamWriter& bam_writer,
	SR_QueryRegion* query_region,
	bam_header_t* bam_header) {

	SR_BamInStreamFree( bam_reader );
	bam_writer.Close();
	SR_QueryRegionFree( query_region );

	// free bam header
	bam_header_destroy( bam_header );
	
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

	// load bam header
	bam_header_t* bam_header = SR_BamInStreamReadHeader( bam_reader );
	
	// bam records are in SR_QueryRegion structure
	SR_QueryRegion* query_region = SR_QueryRegionAlloc();


	while( SR_BamInStreamGetPair( &(query_region->pAnchor), &(query_region->pOrphan), bam_reader ) == SR_OK ) {
		cout << "Got a pair of alignments" << endl;
		bam_writer.SaveAlignment( *query_region->pAnchor );
		bam_writer.SaveAlignment( *query_region->pOrphan );
	}


	// free memory and close files
	Deconstruct( bam_reader, bam_writer, query_region, bam_header );
	cout << "Program done." << endl;

	return 0;

}
