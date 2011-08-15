// This's file is used for prototyping the program,
// and may not be able to be compiled.


int mian( int argc, char** argv ) {

	// ============
	// Declarations
	// ============
	
	// load split-read (SR) candidate to the structure
	// one is a SR anchor record and the other is SR target record
	BamPairRecord bamPairRecord;

	// load reference and corresponding hash table to the structure
	Reference reference;
	HashTable hashTable;

	// query region
	// define the close and far regions
	QueryRegion queryRegion;

	// hash array
	Hashes hashes;

	// interesting type
	SearchModel searchModel;

	// ===================================
	// Load arguments and check parameters
	// ===================================
	
	ParameterParser parameterParser( argc, argv );


	// ============
	// SR alignment
	// ============
	
	while( loadCandidateBamPair( bamPairRecord ) ) { // the function will load candidate ans store it in bamPairRecord
		
		// Are the current reference and hash table correct?
		if ( reference.index < bamPairRecord.anchor.index ) {
			// given the target index, bamPairRecord.index, 
			// load the corresponding reference and hash table in
			jumpReference( reference, bamPairRecord.anchor.index );
			jumpHashTable( hashTable, bamPairRecord.anchor.index );
		} else if ( reference.index > bamPairRecord.index ) {
			// jump bam to the target index, reference.index
			jumpBamFile( reference.index );
		}

		bool splitReadAlignmentFound = false;
		bool stopFinding = false;

		while( !splitReadAlignmentFound && !stopFinding ) {

			// define the search model and store it in searchModel
			setSearchModel( searchModel );
		
			string seqAction;
			if ( searchModel.targetSequenceModel == reverseComplement )
				seqAction = "rc";
			else if ( searchModel.targetSequenceModel == inverse )
				seqAction = "i";
			else
				seqAction = "f";

			// set our target seq according to seqAction
			//*** BY: JIANTAO ***//
			setTargetSequence( bamPairRecord.targetSeq, seqAction.c_str() ); // put our interesting seq in targetSeq

			// set our query region
			//*** BY: JIANTAO ***//
			setQueryRegion( queryRegion, bamPairRecord, fragmentLength );
	
			// get hashes according to query region
			// return an array of the resulting hashes
			//*** BY: JIANTAO ***//
			getHashes( hashes, queryRegion, bamPairRecord.targetSeq );

			// select the hashes that we're interested in
			//=== BY:WAN-PING ===//
			vector< unsigned short > hashIndices;
			selectBestHashes( hashIndex, hashes ); // put the index of our interested hashes in hashIndices
	
			// get alignments for our interested hashes
			//=== BY:WAN-PING ===//
			vector< Alignment > alignments;
			SmithWaterman( alignments, hashIndices, hashes );

			// find SR alignment
			//=== BY:WAN-PING ===//
			vector< Alignment > srAlignments;
			splitReadAlignmentFound = findSplitReadAlignment( srAlignment, alignments, bamPairRecord.target );

			// write the resting SR alignment in bam
			//=== BY:WAN-PING ===//
			if ( splitReadAlignmentFound ) {
				for ( vector< Alignment >::const_iterator ite = srAlignments.begin(); ite != srAlignments.end(); ++ite )
					bamWriter.saveAlignment( *ite );
			} else {
				stopFinding = tryAllPossile();
			}
	}
	



	// ========
	// Clear up
	// ========
	
	free_BamPairRecord( bamPairRecord );
	free_Reference( reference );
	free_HashTable( hashTable );

	bamReader.close();
	bamWriter.close();


	return 0;
}
