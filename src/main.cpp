
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

//#include "utilities/bam_writer.h"
#include "utilities/bam_utilities.h"
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
		bamFile         bam_writer,
		SR_QueryRegion* query_region,
		bam_header_t*   bam_header) {

	SR_BamInStreamFree( bam_reader );
	bam_close( bam_writer );
	SR_QueryRegionFree( query_region );

	// free bam header
	bam_header_destroy( bam_header );
	
}


void CheckFileOrDie( 
		const bamFile bam_writer){

	bool error_found = false;

	if ( bam_writer == NULL ) {
		cout << "ERROR: Cannot open bam file for writing." << endl
		     << "       Please check -o option." << endl;
		error_found = true;
	}

	if ( error_found )
		exit(1);

}

void ResetHeader( bam_header_t* const bam_header ){

	// Reset header line to "@HD\tVN:1.0\tSO:unsorted"
	//BamUtilities::ResetHeaderLineText( bam_header, "@HD\tVN:1.0\tSO:unsorted" );
	// Replace SO:coordinate or SO:queryname by SO:unsorted
	BamUtilities::ReplaceHeaderSoText( bam_header );
}

int main ( int argc, char** argv ) {

	// Parse the arguments and store them
	// The program will exit(1) with printing error message 
	//     if any errors or missing required parameters are found
	const ParameterParser parameter_parser( argc, argv );
	
	// Initialize bam input reader
	// The program will be terminated with printing error message
	//     if the input file cannot be opened.
	SR_BamInStream* bam_reader = SR_BamInStreamAlloc( parameter_parser.input_bam.c_str() );
	// Initialize bam output writer
	bamFile bam_writer = bam_open( parameter_parser.output_bam.c_str(), "w" );

	// Check files statuses
	CheckFileOrDie( bam_writer );

	// Load bam header
	bam_header_t* bam_header = SR_BamInStreamReadHeader( bam_reader );
	if( !parameter_parser.is_input_sorted && !BamUtilities::IsFileSorted( bam_header ) ) {
		// The input bam is unsorted, exit
		cout << "ERROR: The input bam seems unsorted. Please use bamtools sort to sort the bam" << endl
		     << "       or type -s to ignore this checker." << endl;
	}
	
	// Write bam header
	ResetHeader( bam_header );
	bam_header_write( bam_writer, bam_header );

	// bam records are in SR_QueryRegion structure
	SR_QueryRegion* query_region = SR_QueryRegionAlloc();


	for ( unsigned int i = 0; i < 84; ++i ) {
	while( SR_BamInStreamGetPair( &(query_region->pAnchor), &(query_region->pOrphan), bam_reader ) == SR_OK  ) {
		cout << "Got a pair of alignments" << endl;
		bam_write1( bam_writer, query_region->pAnchor );
		bam_write1( bam_writer, query_region->pOrphan );
	}
	}


	// free memory and close files
	Deconstruct( bam_reader, bam_writer, query_region, bam_header );

	cout << "Program done." << endl;

	return 0;

}
