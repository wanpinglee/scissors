#include <string>
#include <stdlib.h>
//#include <stdio.h>
#include <stdint.h>
#include <iostream>

#include "../dataStructures/bam_alignment.h"

using std::string;
using std::cout;
using std::endl;


namespace BamUtilities {

// encode the aligned sequence in bam-foramt encoded sequence
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
				cout << "ERROR: Only the following bases are supported in the BAM format: {=, A, C, G, T, N}." << endl
				     << "       Found: " <<  *pSequence << endl;
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


// Get bam-format packed cigar
bool GetPackedCigar( string& packed_cigar,
		const string& reference,
		const string& query,
		const uint32_t& query_begin,
        	const uint32_t& query_end,
		const uint32_t& read_length) {

	namespace Constant = BamAlignmentConstant;

	// NOTE: the lengths of reference and query should be the same
	uint32_t sequence_length = reference.size();

	if ( query_begin > 0 )
		packed_cigar += query_begin << Constant::kBamCigarShift 
		              | Constant::kBamCsoftClip;

	for ( unsigned int i = 0; i < sequence_length; ++i ) {
		static uint32_t operation_length = 0;
		static uint32_t test_position = i;

		// matchs
		if ( ( reference[i] != '-' ) && ( query[i] != '-' ) ) {
			bool still_go = true;
			while ( still_go ) {
				++operation_length;
				++test_position;
				still_go = ( ( reference[test_position] != '-' ) 
				         && ( query[test_position] != '-' ) 
				         && ( test_position < sequence_length ) );
			}

			packed_cigar += operation_length << Constant::kBamCigarShift 
			              | Constant::kBamCmatch;
		}
		// insertion
		else if ( reference[i] == '-' ) {
			bool still_go = true;
			while( still_go ) {
				++operation_length;
				++test_position;
				still_go = ( ( reference[test_position] == '-' )
				         && ( test_position < sequence_length ) );
			}

			packed_cigar += operation_length << Constant::kBamCigarShift 
			              | Constant::kBamCins; 
		}
		// deletion
		else if ( query[i] == '-' ) {
			bool still_go = true;
			while( still_go ) {
				++operation_length;
				++test_position;
				still_go = ( ( query[i] == '-')
				         && ( test_position < sequence_length ) );
			}

			packed_cigar += operation_length << Constant::kBamCigarShift 
			              | Constant::kBamCdel;
		}
		else {
			cout << "ERROR: Generating the packed cigar failed." << endl;
			cout << "       Unknown char(" << reference[i] << "," << query[i] 
			     << ")." << endl;

			return false;
		}

		i = test_position - 1;

	}

	uint32_t num_tail_soft_clip = read_length - query_end - 1;
	if ( num_tail_soft_clip > 0 )
		packed_cigar += num_tail_soft_clip << Constant::kBamCigarShift
		              | Constant::kBamCsoftClip;


	return true;

}

} // namespace BamUtilities
