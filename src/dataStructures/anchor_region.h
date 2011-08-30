#ifndef DATASTRUCTURES_ANCHOR_REGION_H_
#define DATASTRUCTURES_ANCHOR_REGION_H_

#include <stdint.h>

class AnchorRegion {
 public:
  AnchorRegion() : begin_(0), end_(0){};

  bool IsNewRegion(const uint32_t* packed_cigar, const uint32_t cigar_length, const uint32_t begin);

 private:
  uint32_t CalculateEndPosition(const uint32_t* packed_cigar, const uint32_t cigar_length, const uint32_t begin);
  uint32_t begin_;
  uint32_t end_;
};

#endif // DATASTRUCTURES_ANCHOR_REGION_H_
