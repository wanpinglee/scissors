#include "anchor_region.h"

#include <iostream>

#include "utilities/bam/bam_alignment.h"

const static uint32_t kBamCigarMask = 0x000f;

uint32_t AnchorRegion::CalculateEndPosition(const uint32_t* packed_cigar, const uint32_t packed_cigar_length,
    const uint32_t begin) {
  if (packed_cigar_length == 0)
    return begin;
  
  namespace Constant = BamCigarConstant;
  uint32_t end = begin;
  
  for (uint32_t i = 0; i < packed_cigar_length; ++i) {
    unsigned char op = packed_cigar[i] & kBamCigarMask;
    
    if (op == Constant::kBamCmatch || op == Constant::kBamCdel || op == Constant::kBamCrefSkip)
    	end += packed_cigar[i] >> Constant::kBamCigarShift;
  }
 
  --end;
  return end;
}

bool AnchorRegion::IsNewRegion(const uint32_t* packed_cigar, const uint32_t packed_cigar_length,
    const uint32_t begin) {

  uint32_t end = CalculateEndPosition(packed_cigar, packed_cigar_length, begin);

  // If the current region is a point, we don't consider overlap.
  if (begin_ == end_) {
  	begin_ = begin;
	end_ = end;
	return true;
  }
  
  bool is_begin_overlapped = (begin >= begin_) && (begin <= end_);

  if (is_begin_overlapped) {
    if (end > end_) end_ = end; // assign larger end position
    return false;
  
  } else {
    begin_ = begin;
    end_ = end;
    return true;
  }
}
