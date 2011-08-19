
#include <stdlib.h>

#include <iostream>
#include <string>

extern "C" {
#include "hasher/reader/SR_InHashTable.h"
#include "hasher/reader/SR_Reference.h"
#include "hasher/reader/SR_HashRegionTable.h"
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
	
	// bam file writer
	BamWriter bam_writer( parameter_parser.output_bam );
	bam_writer.Open();


/*
	// check all the I/O
	// By Wan-Ping
	if ( !checkFiles( parameterparser ) ) {
		cout << "ERROR: Errors in I/O are found." << endl;
		return 1;
	}

	// By Jiantao
	// store an alignment
	BamAlignment anchorAlignment, targetAlignment;
	
	// By Jiantao
	// this checker can be put in checkFiles
	BamReader bamReader;
	bamReader.open( parameterparser.inputBam.c_str() );
	if ( !bamReader.isGood() ) {
		cout << "ERROR: Cannot open input bam" << endl;
		return 1;
	}

	// By Wan-Ping
	// this checker can be put in checkFile
	BamWriter bamWriter;
	bamWriter.open( parameterparser.outputBam.c_str() );
	if ( !bamWriter.isGood() ) {
		cout << "ERROR: Cannot open output bam" << endl;
		return 1;
	}

	Reference reference;
	HashTable hashTable;
	loadReference( reference );
	loadHashTable( hashTable );

	while ( loadPairAlignments( anchorAlignment, targetAlignment) ) {  // By Jiantao
	
		if ( reference.index < anchorAlignment.index ) {
			loadReference( reference );
			loadHashTable( hashTable );
		}

		// By Jiantao
		HashRegion* hashRegions;
		// get hash regions for targetAlignment
		getHashRegion( reference, hashTable, anchorAlignment.referenceBegin, targetAlignment.bases, hashRegions );

		// By Wan-Ping
		vector< Alignment > alignments;
		// interface between C and C++
		pickHashRegion( hashRegions, alignments ); // pick hash regions and store them in alignments

		// By Jiantao
		free_HashRegion( hashRegions );

		// By Wan-Ping
		smithWatermanAlign( alignments, reference );

		// By Wan-Ping
		vecotr< Alignment > splitReadAlignments;
		comfirmSplitAlignment( anchorAlignment, alignments, splitReadAlignments ); // refine alignments and decide the result
		BamWriter( anchorAlignment, splitReadAlignments );
	}


	// free memory
	// By Jiantao
	free_BamAlignment( anchorAlignment );
	free_BamAlignment( targetAlignment );
	free_Reference( reference );
	free_HashTable( hashTable );
	bamReader.close();
	// By Wan-Ping
	bamWriter.close();

	
	// ================
	// END of prototype
	// ================
*/


/*
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
*/

	bam_writer.Close();
	return 0;

}
