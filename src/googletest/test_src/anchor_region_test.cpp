#include "dataStructures/anchor_region.h"

#include <stdint.h>

#include "gtest/gtest.h"

#include "utilities/bam/bam_alignment.h"

namespace {
// convert the readable cigar string into the packed string
void ConvertCigarToPackedCigar( vector<uint32_t>* packed_cigar, const string& cigar){
	namespace Constant = BamAlignmentConstant;
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
				printf("ERROR: A cigar operator should follow a number.\n");
			else {
				current_cigar = atoi( cigar.substr( digit_begin, digit_length ).c_str() ) << Constant::kBamCigarShift;
				
				switch ( cigar[i] ) {
					case 'S': current_cigar |= Constant::kBamCsoftClip; break;
					case 'I': current_cigar |= Constant::kBamCins;      break;
					case 'D': current_cigar |= Constant::kBamCdel;      break;
					case 'M': current_cigar |= Constant::kBamCmatch;    break;
					default:
					  printf("ERROR: An unknown cigar operator, %c.\n", cigar[i]);
				}
			}
			
			packed_cigar->push_back( current_cigar );
			
			digit_found  = false;
			digit_length = 0;
		}
	}
}

bool ConvertFromVectorToArray(const vector<uint32_t>& i_vector, uint32_t*& o_array, uint32_t* array_size ) {
  // clean o_array
  if ((o_array != NULL) && (*array_size > 1))
    delete [] o_array;
  else if ((o_array != NULL) && (*array_size == 1))
    delete o_array;

  *array_size = i_vector.size();
  o_array = new uint32_t [i_vector.size()];
  for (uint32_t i = 0; i < *array_size; ++i) {
    o_array[i] = i_vector[i];
  }

  return true;
}

TEST(AnchorRegionTest, Test) {
  AnchorRegion anchor_region;

  vector<uint32_t> packed_cigar;
  uint32_t* packed_cigar_array = NULL;
  uint32_t  packed_cigar_array_size = 0;

  bool new_region;
  // first anchor region
  ConvertCigarToPackedCigar(&packed_cigar, "100M");
  ConvertFromVectorToArray(packed_cigar, packed_cigar_array, &packed_cigar_array_size);
  new_region = anchor_region.IsNewRegion(packed_cigar_array, packed_cigar_array_size, 0);
  ASSERT_TRUE(new_region);

  // second anchor region
  ConvertCigarToPackedCigar(&packed_cigar, "90M");
  ConvertFromVectorToArray(packed_cigar, packed_cigar_array, &packed_cigar_array_size);
  new_region = anchor_region.IsNewRegion(packed_cigar_array, packed_cigar_array_size, 0);
  ASSERT_FALSE(new_region);

  ConvertCigarToPackedCigar(&packed_cigar, "20M5D60M");
  ConvertFromVectorToArray(packed_cigar, packed_cigar_array, &packed_cigar_array_size);
  new_region = anchor_region.IsNewRegion(packed_cigar_array, packed_cigar_array_size, 1000);
  ASSERT_TRUE(new_region);

  ConvertCigarToPackedCigar(&packed_cigar, "100M");
  ConvertFromVectorToArray(packed_cigar, packed_cigar_array, &packed_cigar_array_size);
  new_region = anchor_region.IsNewRegion(packed_cigar_array, packed_cigar_array_size, 1050);
  ASSERT_FALSE(new_region);

  ConvertCigarToPackedCigar(&packed_cigar, "100M");
  ConvertFromVectorToArray(packed_cigar, packed_cigar_array, &packed_cigar_array_size);
  new_region = anchor_region.IsNewRegion(packed_cigar_array, packed_cigar_array_size, 1150);
  ASSERT_TRUE(new_region);


  // free packed_cigar_array
  if ((packed_cigar_array != NULL) && (packed_cigar_array_size > 1))
    delete [] packed_cigar_array;
  else if ((packed_cigar_array != NULL) && (packed_cigar_array_size == 1))
    delete packed_cigar_array;
}

} // namespace
