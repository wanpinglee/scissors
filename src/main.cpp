
#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
//#include "hasher/reader/SR_InHashTable.h"
#include "hasher/reader/SR_Reference.h"
//#include "hasher/reader/SR_HashRegionTable.h"
#include "utilities/SR_BamInStream.h"
#include "hasher/common/SR_Types.h"
#include "dataStructures/SR_QueryRegion.h"
#include "samtools/bam.h"
}

//#include "utilities/bam_writer.h"
#include "utilities/bam_utilities.h"
#include "utilities/parameter_parser.h"

using std::string;
using std::cout;
using std::endl;

struct MainFiles {
	SR_BamInStream* bam_reader;  // bam reader
	bamFile         bam_writer;  // bam writer
	FILE*           ref_reader;  // reference reader
	FILE*           hash_reader; // hash table reader
};

struct MainVars{
	SR_QueryRegion* query_region;
	bam_header_t*   bam_header;
	SR_Reference*   reference;
};


void Deconstruct( MainFiles& files, MainVars& vars ) {
	// close files
	SR_BamInStreamFree( files.bam_reader );
	bam_close( files.bam_writer );
	fclose( files.ref_reader );
	fclose( files.hash_reader );

	// free variables
	SR_QueryRegionFree( vars.query_region );
	bam_header_destroy( vars.bam_header );
	SR_ReferenceFree( vars.reference );

	
}

void CheckFileOrDie( 
		const ParameterParser& parameter_parser,
		const MainFiles& files){

	bool error_found = false;

	if ( files.bam_writer == NULL ) {
		cout << "ERROR: Cannot open " << parameter_parser.output_bam << " for writing." << endl
		     << "       Please check -o option." << endl;
		error_found = true;
	}

	if ( files.ref_reader == NULL ) {
		cout << "ERROR: Cannot open " << parameter_parser.reference_filename << " for reading." << endl
		     << "       Please check -r option." << endl;
		error_found = true;
	}

	if ( files.hash_reader == NULL ) {
		cout << "ERROR: Cannot open " << parameter_parser.hash_filename << " for reading." << endl
		     << "       Please check -r option." << endl;
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


	// =================
	// Files Preparation
	// =================
	MainFiles files;
	
	// Initialize bam input reader
	// The program will be terminated with printing error message
	//     if the input file cannot be opened.
	files.bam_reader = SR_BamInStreamAlloc( parameter_parser.input_bam.c_str() );
	// Initialize bam output writer
	files.bam_writer = bam_open( parameter_parser.output_bam.c_str(), "w" );

	// Initialize reference input reader
	files.ref_reader = fopen( parameter_parser.reference_filename.c_str(), "rb");
	files.hash_reader = fopen( parameter_parser.hash_filename.c_str(), "rb");

	// Check files statuses
	CheckFileOrDie( parameter_parser, files );


	// =====================
	// Variables Preparation
	// =====================
	MainVars vars;

	// Load bam header
	vars.bam_header = SR_BamInStreamReadHeader( files.bam_reader );
	if( !parameter_parser.is_input_sorted && !BamUtilities::IsFileSorted( vars.bam_header ) ) {
		// The input bam is unsorted, exit
		cout << "ERROR: The input bam seems unsorted. "
		     << "Please use bamtools sort to sort the bam" << endl
		     << "       or type -s to ignore this checker." << endl;
		exit(1);
	}
	
	// Write bam header
	ResetHeader( vars.bam_header );
	bam_header_write( files.bam_writer, vars.bam_header );


	// =====================
	// Variables Preparation
	// =====================
	
	// bam records are in SR_QueryRegion structure
	vars.query_region = SR_QueryRegionAlloc();

	// Load reference
	uint32_t buffer_size = 1000000;
	vars.reference = SR_ReferenceAlloc( buffer_size );
	SR_ReferenceRead( vars.reference, files.ref_reader );

	// =========
	// Algorithm
	// =========
	
	for ( unsigned int i = 0; i < 84; ++i ) {
	while( SR_BamInStreamGetPair( &(vars.query_region->pAnchor), &(vars.query_region->pOrphan), files.bam_reader ) == SR_OK  ) {
		cout << "Got a pair of alignments" << endl;
		bam_write1( files.bam_writer, vars.query_region->pAnchor );
		bam_write1( files.bam_writer, vars.query_region->pOrphan );
	}
	}


	// free memory and close files
	Deconstruct( files, vars );

	cout << "Program done." << endl;

	return 0;

}
