#include "anchor_region.h"

#include "dataStructures/bam_alignment.h"

const static unsigned char kBamCigarMask = 0x0f;

uint32_t AnchorRegion::CalculateEndPosition(const char* cigar, const uint32_t cigar_length,
    const uint32_t begin) {
  namespace Constant = BamAlignmentConstant;
  uint32_t end = begin;
  
  for (uint32_t i = 0; i < cigar_length; ++i) {
    unsigned char op = cigar[i] & kBamCigarMask;
    if (op == Constant::kBamCmatch || op == Constant::kBamCdel || op == Constant::kBamCrefSkip)
    end += cigar[i] >> Constant::kBamCigarShift;
  }
  
  return end;
}

bool AnchorRegion::IsNewRegion(const char* cigar, const uint32_t cigar_length,
    const uint32_t begin) {

  bool is_begin_overlapped = (begin >= begin_) && (begin <= end_);
  if (is_begin_overlapped) {
    uint32_t end = CalculateEndPosition(cigar, cigar_length, begin);
    if (end > end_) end_ = end; // assign larger end position
    return true;
  
  } else {
    return false;
  }
}
