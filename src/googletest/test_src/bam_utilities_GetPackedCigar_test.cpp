#include <iostream>
#include <vector>

#include "gtest/gtest.h"

#include "utilities/bam/bam_utilities.h"
#include "utilities/bam/bam_constant.h"

using std::cout;
using std::endl;
using std::vector;

namespace {

// convert the readable cigar string into the packed string
void ConvertCigarToPackedCigar( vector<uint32_t>* packed_cigar, const string& cigar){
	namespace Constant = BamCigarConstant;
	packed_cigar->clear();
	uint32_t digit_begin  = 0;
	uint32_t digit_length = 0;
	bool digit_found = false;

	for ( unsigned int i = 0; i < cigar.size(); ++i ) {
		if ( isdigit( cigar[i] ) ) {
			if ( !digit_found )
				digit_begin = i;
			
			++digit_length;
			digit_found = true;
		} else {
			uint32_t current_cigar = 0;
			if ( !digit_found )
				cout << "ERROR: A cigar operator should follow a number." << endl;
			else {
				current_cigar = atoi( cigar.substr( digit_begin, digit_length ).c_str() ) << Constant::kBamCigarShift;
				
				switch ( cigar[i] ) {
					case 'S': current_cigar |= Constant::kBamCsoftClip; break;
					case 'I': current_cigar |= Constant::kBamCins;      break;
					case 'D': current_cigar |= Constant::kBamCdel;      break;
					case 'M': current_cigar |= Constant::kBamCmatch;    break;
					default:
					  cout << "ERROR: An unknown cigar operator, " << cigar[i] << "." << endl;
				}
			}
			
			packed_cigar->push_back( current_cigar );
			
			digit_found  = false;
			digit_length = 0;
		}
	}
}

TEST(GetPackedCigarTest, NullString) {
  string reference, query;
  uint32_t query_begin = 0;
  uint32_t query_end   = 0;
  uint32_t read_length = 50;
  vector<uint32_t> packed_cigar;

  // ===============
  // expected result
  // ===============
  vector<uint32_t> expected_cigar;
  ConvertCigarToPackedCigar( &expected_cigar, "50S");
  // ===============

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
  	query_begin, query_end, read_length);


  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );

}

TEST(GetPackedCigarTest, Match) {
  // ===========
  // Test case 1
  // ===========
  string reference = "ACTGGTACGTGTCAGTCGTACCAGTCTGAC";
  string query     = "ACTGGTACGTGTCAGTCGTACCAGTCTGAC";
  uint32_t query_begin = 0;
  uint32_t query_end   = 29;
  uint32_t read_length = 30;
  vector<uint32_t> packed_cigar;
  

  vector<uint32_t> expected_cigar;
  ConvertCigarToPackedCigar( &expected_cigar, "30M");
	
  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
	query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


  // ===========
  // Test case 2
  // ===========
  reference = "ACTGGTACGTGTCAGTCGTACCAGTCTGAC";
  query     = "ACTGGTACGTGTCAGTCGTACCAGTCTGAC";
  query_begin = 0;
  query_end   = 29;
  read_length = 50;

  ConvertCigarToPackedCigar( &expected_cigar, "30M20S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
	query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


  // ===========
  // Test case 3
  // ===========
  reference = "ACTGGTACGTGTCAGTCGTACCAGTCTGAC";
  query     = "GGGGGTACGTGTCAGTCGTACCAGTCTTTT";
  query_begin = 10;
  query_end   = 39;
  read_length = 50;

  ConvertCigarToPackedCigar( &expected_cigar, "10S30M10S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
	query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );

  // ===========
  // Test case 4
  // ===========
  reference = "ACTGGTACGTGTCAGTC";
  query     = "GGGGGTACGTGTCAGTC";
  query_begin = 0;
  query_end   = 16;
  read_length = 40;

  ConvertCigarToPackedCigar( &expected_cigar, "17M23S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
	query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i ) 
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


}

TEST(GetPackedCigarTest, OneBaseMatch) {
  // ===========
  // Test case 1
  // ===========
  string reference = "A";
  string query     = "A";
  uint32_t query_begin = 0;
  uint32_t query_end   = 0;
  uint32_t read_length = 1;
  vector<uint32_t> packed_cigar;

  // expected result
  vector<uint32_t> expected_cigar;
  ConvertCigarToPackedCigar( &expected_cigar, "1M");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


  // ===========
  // Test case 2
  // ===========
  reference = "A";
  query     = "A";
  query_begin = 0;
  query_end   = 0;
  read_length = 3;

  // expected result
  ConvertCigarToPackedCigar( &expected_cigar, "1M2S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


  // ===========
  // Test case 3
  // ===========
  reference = "T";
  query     = "T";
  query_begin = 2;
  query_end   = 2;
  read_length = 3;

  // expected result
  ConvertCigarToPackedCigar( &expected_cigar, "2S1M");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


  // ===========
  // Test case 4
  // ===========

  reference = "C";
  query     = "C";
  query_begin = 1;
  query_end   = 1;
  read_length = 3;

  // expected result
  ConvertCigarToPackedCigar( &expected_cigar, "1S1M1S");


  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( packed_cigar.size(), expected_cigar.size() );

}

TEST(GetPackedCigarTest, Insertion) {
  // ===========
  // Test case 1
  // ===========
  string reference = "AC-TTG";
  string query     = "ACGTTG";
  uint32_t query_begin = 0;
  uint32_t query_end   = 5;
  uint32_t read_length = 6;
  vector<uint32_t> packed_cigar;

  // expected result
  vector<uint32_t> expected_cigar;
  ConvertCigarToPackedCigar( &expected_cigar, "2M1I3M");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


  // ===========
  // Test case 2
  // ===========
  reference = "AC-T-G";
  query     = "ACGTTG";
  query_begin = 0;
  query_end   = 5;
  read_length = 6;

  // expected result
  ConvertCigarToPackedCigar( &expected_cigar, "2M1I1M1I1M");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );

  // ===========
  // Test case 3
  // ===========
  reference = "-C-T-G";
  query     = "ACGTTG";
  query_begin = 3;
  query_end   = 8;
  read_length = 13;

  // expected result
  ConvertCigarToPackedCigar( &expected_cigar, "3S1I1M1I1M1I1M4S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );
}

TEST(GetPackedCigarTest, Deletion) {
  // ===========
  // Test case 1
  // ===========
  string reference = "ACGTTG";
  string query     = "AC-TTG";
  uint32_t query_begin = 0;
  uint32_t query_end   = 4;
  uint32_t read_length = 6;
  vector<uint32_t> packed_cigar;

  // expected result
  vector<uint32_t> expected_cigar;
  ConvertCigarToPackedCigar( &expected_cigar, "2M1D3M1S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );


  // ===========
  // Test case 2
  // ===========
  reference = "ACGTTG";
  query     = "AC---G";
  query_begin = 2;
  query_end   = 4;
  read_length = 10;

  // expected result
  ConvertCigarToPackedCigar( &expected_cigar, "2S2M3D1M5S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );

}

TEST(GetPackedCigarTest, Mix) {
  // ===========
  // Test case 1
  // ===========
  string reference = "ACGT--GATA";
  string query     = "A--TTAGA-A";
  uint32_t query_begin = 3;
  uint32_t query_end   = 9;
  uint32_t read_length = 12;
  vector<uint32_t> packed_cigar;

  // expected result
  vector<uint32_t> expected_cigar;
  ConvertCigarToPackedCigar( &expected_cigar, "3S1M2D1M2I2M1D1M2S");

  Scissors::BamUtilities::GetPackedCigar( packed_cigar, reference, query,
        query_begin, query_end, read_length);

  for ( unsigned int i = 0; i < expected_cigar.size(); ++i )
  	EXPECT_EQ( expected_cigar[i], packed_cigar[i] );
  EXPECT_EQ( expected_cigar.size(), packed_cigar.size() );
}
} // namespace
