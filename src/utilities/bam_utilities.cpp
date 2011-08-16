#include <string>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


void EncodeQuerySequence( string& encodedSequence, const string& sequence ) {
	
	// prepare the encoded query string
	const unsigned int sequenceLen = sequence.size();
	const unsigned int encodedSequenceLen = (unsigned int)( ( sequenceLen / 2.0 ) + 0.5 );
	encodedSequence.resize( encodedSequenceLen );
	char* pEncodedSequence = (char*) encodedSequence.data();
	const char* pSequence = (const char*) sequence.data();

	unsigned char nucleotideCode;
	bool useHighWord = true;

	for( unsigned int i = 0; i < sequenceLen; ++i ) {
		switch( *pSequence ) {
			case '=':
				nucleotideCode = 0;
				break;
			case 'A':
				nucleotideCode = 1;
				break;
			case 'C':
				nucleotideCode = 2;
				break;
			case 'G':
				nucleotideCode = 4;
				break;
			case 'T':
				nucleotideCode = 8;
				break;
			case 'N':
				nucleotideCode = 15;
				break;
			default:
				printf( "ERROR: Only the following bases are supported in the BAM format: {=, A, C, G, T, N}. Found [%c]\n", *pSequence );
				exit(1);
		}

		// pack the nucleotide code
		if( useHighWord ) {
			*pEncodedSequence = nucleotideCode << 4;
			useHighWord = false;
		} else {
			*pEncodedSequence |= nucleotideCode;
			++pEncodedSequence;
			useHighWord = true;
		}

		// increment the query position
		++pSequence;
	}
}
