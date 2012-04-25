#include "bam_utilities.h"

#include <string>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <sstream>

#include "outsources/samtools/sam_header.h"
#include "bam_constant.h"

using std::string;
using std::cerr;
using std::endl;

static const char COMPLEMENT_BASE_MAP[16] = {0x0, 0x8, 0x4, 0xf, 0x2,0xf,0xf,0xf,0x1,0xf,0xf,0xf,0xf,0xf,0xf,0xf};

namespace Scissors {
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

bool IsFileSorted( const bam_header_t* const header ){
	string header_string( header->text );
	size_t so_pos = header_string.find("SO:coordinate");

	if ( so_pos != string::npos )
		return true;
	else
		return false;
}

bool ReplaceHeaderSoText(bam_header_t* const header) {
	
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

	namespace Constant = BamCigarConstant;

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
			cerr << "ERROR: Generating the packed cigar failed." << endl;
			cerr << "       Unknown char(" << reference[i] << "," << query[i] 
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

void GetReverseComplementSequence(const uint8_t* original, 
                                  const int& length, 
              	 		  uint8_t* reverse) 
{
  int length_even = ((length % 2) == 1) ? length - 1 : length;
  for (int i = 0; i < length_even; i=i+2) {
    char upper = COMPLEMENT_BASE_MAP[bam1_seqi(original, length - i - 1)];
    char lower = COMPLEMENT_BASE_MAP[bam1_seqi(original, length - i - 2)];
    *reverse = (upper << 4) | lower;
    ++reverse;
  }

  if ((length % 2) == 1) {
    char upper = COMPLEMENT_BASE_MAP[bam1_seqi(original, 0)];
    *reverse = (upper << 4) & 0xf0;
    ++reverse;
  }
}

void GetComplementSequence(const uint8_t* original,
                           const int& length,
			   uint8_t* reverse)
{
  for (int i = 0; i < length; i=i+2) {
    char upper = COMPLEMENT_BASE_MAP[bam1_seqi(original, i)];
    char lower = COMPLEMENT_BASE_MAP[bam1_seqi(original, i+1)];
    *reverse = (upper << 4) | lower;
    ++reverse;
  }
}

void GetInverseSequence(const uint8_t* original,
                        const int& length,
			uint8_t* reverse)
{
  int length_even = ((length % 2) == 1) ? length - 1 : length;
  for (int i = 0; i < length_even; i=i+2) {
    char upper = bam1_seqi(original, length - i - 1);
    char lower = bam1_seqi(original, length - i - 2);
    *reverse = (upper << 4) | lower;
    ++reverse;
  }

  if ((length % 2) == 1) {
    char upper = bam1_seqi(original, 0);
    *reverse = (upper << 4) & 0xf0;
    ++reverse;
  }
}

// TODO @ Wan-Ping: isize, flag, and mq are not yet to be assigned.
void ConvertAlignmentToBam1(const Alignment& al, 
                            const bam1_t& original_record, 
			    bam1_t* new_record) {
  //cout << original_record.core.l_qname << endl;
  //cout << al.query_begin << "\t" << al.query_end << endl;
  vector<uint32_t> packed_cigar;
  GetPackedCigar(packed_cigar, 
                 al.reference, 
		 al.query, 
		 al.query_begin, 
		 al.query_end, 
		 original_record.core.l_qseq);

  // copy the core info
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

  // set the flag
  if (al.is_seq_inverse && al.is_seq_complement)
    new_record->core.flag |= 0x0010;

  // so far, no extra aux info is allowed
  new_record->l_aux = 0;

  // calculate length of data
  int data_length = new_record->core.l_qname +
                    new_record->core.n_cigar * 4 +
		    (new_record->core.l_qseq + 1) / 2 +
		    new_record->core.l_qseq +
		    new_record->l_aux;

  // set data length
  new_record->data_len = data_length;
  new_record->m_data   = data_length;
  kroundup32(new_record->m_data);

  // set data
  uint8_t* data = (uint8_t*) calloc(new_record->m_data, sizeof(uint8_t));  // Thread.cpp will delete those
  uint8_t* data_ptr = data;
  // copy the read name
  memcpy(data_ptr, original_record.data, new_record->core.l_qname);
  data_ptr += new_record->core.l_qname;

  // set the cigar
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

  if (!al.is_seq_inverse && !al.is_seq_complement) { // forward
    memcpy(data_ptr, bam1_seq(&original_record), (new_record->core.l_qseq + 1) / 2);
  } else if (al.is_seq_inverse && al.is_seq_complement) { // reverse complement
    GetReverseComplementSequence(bam1_seq(&original_record), new_record->core.l_qseq, data_ptr);
  } else if (!al.is_seq_inverse && al.is_seq_complement) { // complement
    GetComplementSequence(bam1_seq(&original_record), new_record->core.l_qseq, data_ptr);
  } else { // inverse
    GetInverseSequence(bam1_seq(&original_record), new_record->core.l_qseq, data_ptr);
  }

  data_ptr += ((new_record->core.l_qseq + 1) / 2);

  // copy base qualities
  memcpy(data_ptr, bam1_qual(&original_record), new_record->core.l_qseq);
  data_ptr += new_record->core.l_qseq;

  new_record->data = data;


}

bool AppendReferenceSequence(const char** names,
                             const uint32_t* lens,
			     const char** md5s,
                             const int& n_sequences,
			     bam_header_t* const header) {
  // for sizing
  char char_a;
  uint32_t uint32_t_a;

  int total_names = header->n_targets + n_sequences;
  char** new_names = (char**)calloc(total_names, sizeof(&char_a));
  uint32_t* new_lens = (uint32_t*)calloc(total_names, sizeof(uint32_t_a));

  // copy the original names and lengths
  for (int i = 0; i < header->n_targets; ++i) {
    int len_name = strlen(header->target_name[i]);
    new_names[i] = (char*)calloc(len_name, sizeof(char_a));
    memcpy(new_names[i], header->target_name[i], len_name);
    new_lens[i] = header->target_len[i];
    free(header->target_name[i]);
  }
  free(header->target_name);
  free(header->target_len);

  std::ostringstream additional_text;
  additional_text << endl;
  // append the new names and lengths
  for (int i = 0; i < n_sequences; ++i) {
    int len_name = strlen(names[i]);
    new_names[i + header->n_targets] = (char*)calloc(len_name, sizeof(char_a));
    memcpy(new_names[i + header->n_targets], names[i], len_name);
    new_lens[i + header->n_targets] = lens[i];
    additional_text << "@SQ\tSN:" << names[i] << "\tLN:" << lens[i] << "\tM5:" << md5s[i] << endl;
  }
  header->n_targets += n_sequences;
  header->target_name = new_names;
  header->target_len  = new_lens;

  string original_text = header->text;
  size_t sq_begin = original_text.find("@SQ");
  size_t sq_end   = sq_begin + 1;
  while (true) {
    size_t found = original_text.find("@SQ", sq_end);
    if (found != string::npos) { // found
      sq_end = found + 1;
    } else { // no more @SQ
      found = original_text.find("\n", sq_end);
      sq_end = found;
      break;
    }
  }
  
  string new_text = original_text.substr(0, sq_end);
  new_text += additional_text.str();
  new_text += original_text.substr(sq_end + 1);
  free(header->text);
  header->text = (char*) calloc(new_text.size(), sizeof(char_a));
  memcpy(header->text, new_text.c_str(), new_text.size());
  header->l_text = new_text.size();

  return true;
}

} // namespace BamUtilities
} //namespace Scissors
