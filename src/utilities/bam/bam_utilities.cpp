#include "bam_utilities.h"

#include <string>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <sstream>

#include "outsources/samtools/sam_header.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

namespace BamAlignmentConstant {

const unsigned char kBamCoreSize   = 32;
const unsigned char kBamCigarShift = 4;

// cigar operator for bam packed cigar
const unsigned char kBamCmatch     = 0;
const unsigned char kBamCins       = 1;
const unsigned char kBamCdel       = 2;
const unsigned char kBamCrefSkip   = 3;
const unsigned char kBamCsoftClip  = 4;
const unsigned char kBamChardClip  = 5;
const unsigned char kBamCpad       = 6;

} // namespace BamAlignmentConstant


namespace BamUtilities {
bool ResetHeaderText( bam_header_t* const header, const string& header_string ){
	
	// Follow by the samtools style to allocate memory
	if ( header->text ) {
		header->l_text = header_string.size();
		free( header->text );
		header->text = (char*)calloc(header->l_text + 1, 1);
		strncpy( header->text, header_string.c_str(), header->l_text );
	}

	return true;
}
/*
bool ResetHeaderLineText( bam_header_t* const header, const char* const str ) {
	string input( header->text );

	// The header line (@HD) must be the first line if present
	if ( input.empty() )
		input = str;
	else {
		if ( input.substr(1, 2) != "HD" ) {
		// @HD doesn't exist
			input.insert(0, str+'\n');
		}
		else {
		// @HD is found
			size_t pos = input.find('\n');
			if ( pos == string::npos ) // cannot find '\n'
				pos = input.find('\r');

			if ( pos == string::npos ) // cannot find '\n' and '\r'
				return false;


		}
	}

	ResetHeaderText( header, input );
	return true;
}
*/

bool IsFileSorted( const bam_header_t* const header ){
	string header_string( header->text );
	size_t so_pos = header_string.find("SO:coordinate");

	if ( so_pos != string::npos )
		return true;
	else
		return false;
}

bool ReplaceHeaderSoText( bam_header_t* const header ) {
	
	string header_string( header->text );
	
	size_t so_pos = header_string.find("SO:coordinate");
	if ( so_pos != string::npos )
		header_string.replace(so_pos, 13, "SO:unsorted");

	so_pos = header_string.find("SO:queryname");
	if ( so_pos != string::npos )
		header_string.replace(so_pos, 12, "SO:unsorted");

	ResetHeaderText(header, header_string);

	return true;

}

// encode the aligned sequence in bam-foramt encoded sequence
void EncodeQuerySequence( string& encodedSequence, const string& sequence ) {
	
	// prepare the encoded query string
	const unsigned int sequenceLen = sequence.size();
	const unsigned int encodedSequenceLen = (unsigned int)( ( sequenceLen / 2.0 ) + 0.5 );
	encodedSequence.resize( encodedSequenceLen );
	char* pEncodedSequence = (char*) encodedSequence.c_str();
	const char* pSequence = (const char*) sequence.c_str();

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
				++pSequence;
				continue;
				//cout << "ERROR: Only the following bases are supported in the BAM format: {=, A, C, G, T, N}." << endl
				//     << "       Found: " <<  *pSequence << endl;
				//exit(1);
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
bool GetPackedCigar( vector<uint32_t>& packed_cigar,
		const string& reference,
		const string& query,
		const uint32_t& query_begin,
        	const uint32_t& query_end,
		const uint32_t& read_length) {

	namespace Constant = BamAlignmentConstant;

	packed_cigar.clear();


	// NOTE: the lengths of reference and query should be the same
	uint32_t sequence_length = reference.size();
	uint32_t current_cigar   = 0;
	
	// no aligned bases
	if ( sequence_length == 0 ) {
		current_cigar = read_length << Constant::kBamCigarShift
		              | Constant::kBamCsoftClip;
		packed_cigar.push_back(current_cigar);
		return true;
	}
	

	// beginning soft clips
	if ( query_begin > 0 ) {
		current_cigar = query_begin << Constant::kBamCigarShift 
		              | Constant::kBamCsoftClip;
		packed_cigar.push_back(current_cigar);

	}

	//const unsigned int test_length = query_end - query_begin + 1;
	for ( unsigned int i = 0; i < sequence_length; ++i ) {
		uint32_t operation_length = 0;
		uint32_t test_position = i;

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

			current_cigar = operation_length << Constant::kBamCigarShift 
			              | Constant::kBamCmatch;
		} // if
		// insertion
		else if ( reference[i] == '-' ) {
			bool still_go = true;
			while( still_go ) {
				++operation_length;
				++test_position;
				still_go = ( ( reference[test_position] == '-' )
				         && ( test_position < sequence_length ) );
			}

			current_cigar = operation_length << Constant::kBamCigarShift 
		                      | Constant::kBamCins; 
		} // else if
		// deletion
		else if ( query[i] == '-' ) {
			bool still_go = true;
			while( still_go ) {
				++operation_length;
				++test_position;
				still_go = ( ( query[test_position] == '-')
				         && ( test_position < sequence_length ) );
			}

			current_cigar = operation_length << Constant::kBamCigarShift 
		     	              | Constant::kBamCdel;
		} // else if
		else {
			cout << "ERROR: Generating the packed cigar failed." << endl;
			cout << "       Unknown char(" << reference[i] << "," << query[i] 
			     << ")." << endl;

			return false;
		} // else

		i = test_position - 1;
		packed_cigar.push_back(current_cigar);

	} // for

	// trailing soft clips
	uint32_t num_tail_soft_clip = read_length - query_end - 1;
	if ( num_tail_soft_clip > 0 ) {
		current_cigar = num_tail_soft_clip << Constant::kBamCigarShift
		              | Constant::kBamCsoftClip;
		packed_cigar.push_back(current_cigar);

	}


	return true;

}

// TODO @ Wan-Ping: isize is not yet to be assigned.
void ConvertAlignmentToBam1(const Alignment& al, 
                            const bam1_t& original_record, 
			    bam1_t* new_record) {
  cout << original_record.core.l_qname << endl;
  cout << al.query_begin << "\t" << al.query_end << endl;
  vector<uint32_t> packed_cigar;
  GetPackedCigar(packed_cigar, 
                 al.reference, 
		 al.query, 
		 al.query_begin, 
		 al.query_end, 
		 original_record.core.l_qseq);
  cout << packed_cigar.size() << endl;

  new_record->core.tid     = original_record.core.tid;
  new_record->core.pos     = al.reference_begin;
  new_record->core.bin     = bam_reg2bin(al.query_begin, al.query_end);
  new_record->core.qual    = al.quality;
  new_record->core.l_qname = original_record.core.l_qname;
  new_record->core.flag    = 0;
  new_record->core.n_cigar = packed_cigar.size();
  new_record->core.l_qseq  = original_record.core.l_qseq;
  new_record->core.mtid    = original_record.core.mtid;
  new_record->core.mpos    = original_record.core.mpos;
  new_record->core.isize   = 0;

  new_record->l_aux = 0;

  int data_length = new_record->core.l_qname +
                    new_record->core.n_cigar * 4 +
		    (new_record->core.l_qseq + 1) / 2 +
		    new_record->core.l_qseq +
		    new_record->l_aux;

  new_record->data_len = data_length;
  new_record->m_data   = data_length;

  uint8_t* data = new uint8_t[data_length];  // Thread.cpp will delete those
  uint8_t* data_ptr = data;
  memcpy(data_ptr, original_record.data, new_record->core.l_qname);
  data_ptr += new_record->core.l_qname;

  for (unsigned int i = 0; i < packed_cigar.size(); ++i) {
    *data_ptr = static_cast<uint8_t>(packed_cigar[i]);
    ++data_ptr;
    *data_ptr = static_cast<uint8_t>(packed_cigar[i] >> 8);
    ++data_ptr;
    *data_ptr = static_cast<uint8_t>(packed_cigar[i] >> 16);
    ++data_ptr;
    *data_ptr = static_cast<uint8_t>(packed_cigar[i] >> 24);
    ++data_ptr;
  }

  memcpy(data_ptr, bam1_seq(&original_record), (new_record->core.l_qseq + 1) / 2);
  data_ptr += ((new_record->core.l_qseq + 1) / 2);

  memcpy(data_ptr, bam1_qual(&original_record), new_record->core.l_qseq);
  data_ptr += new_record->core.l_qseq;

  new_record->data = data;


}

} // namespace BamUtilities
