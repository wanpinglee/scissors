#ifndef DATASTRUCTURES_ANCHOR_REGION_H_
#define DATASTRUCTURES_ANCHOR_REGION_H_

#include <stdint.h>

// Short description: 
//   Detect if the current anchor region overlaps with the previous one.
//
// Long description:
//   Each region is described by a begin and an end points; the end points 
//   actually is calculated by the begin point and its bam packed cigar.
//   IsNewRegion is a function to check if the given region is a new region 
//   which means that the given region does *not* overlap with the region 
//   indicated by begin_ and end_. If an overlap is *not* detected, begin_
//   and end_ will be updated to the given region of IsNewRegion. On the 
//   contrary, if an overlap is detected, begin_ and end_ may be updated to 
//   extend its coverage of the given region of IsNewRegion.
//   Note that when begin_ is equal to end_, IsNewRegion always returns true.

namespace Scissors {

class AnchorRegion {
 public:
  AnchorRegion() : begin_(0), end_(0){};

  bool IsNewRegion(const uint32_t* packed_cigar, const uint32_t cigar_length, const uint32_t begin);

 private:
  uint32_t CalculateEndPosition(const uint32_t* packed_cigar, const uint32_t cigar_length, const uint32_t begin);
  uint32_t begin_;
  uint32_t end_;
};
} //namespace Scissors

#endif // DATASTRUCTURES_ANCHOR_REGION_H_
