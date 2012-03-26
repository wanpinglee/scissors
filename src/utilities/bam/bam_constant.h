#ifndef UTILITIES_BAM_BAM_CONSTANT_H_
#define UTILITIES_BAM_BAM_CONSTANT_H_

#include <stdint.h>

namespace BamCigarConstant {
const unsigned char kBamCigarShift = 4;

// cigar operator for bam packed cigar
const unsigned char kBamCmatch     = 0;
const unsigned char kBamCins       = 1;
const unsigned char kBamCdel       = 2;
const unsigned char kBamCrefSkip   = 3;
const unsigned char kBamCsoftClip  = 4;
const unsigned char kBamChardClip  = 5;
const unsigned char kBamCpad       = 6;
} // namespace BamCigarConstant

namespace BamFlagConstant {
const uint16_t kBamFPaired       = 1;
const uint16_t kBamFProperPair   = 2;
const uint16_t kBamFUnmapped     = 4;
const uint16_t kBamFUnmappedMate = 8;
const uint16_t kBamFReverse      = 16;
const uint16_t kBamFReverseMate  = 32;
const uint16_t kBamFMate1        = 64;
const uint16_t kBamFMate2        = 128;
const uint16_t kBamFNotPrimary   = 256;
const uint16_t kBamFQcFail       = 512;
const uint16_t kBamFDuplication  = 1024;
} // namespace BamFlagConstant

namespace BamOtherConstant {
const unsigned char kBamCoreSize   = 32;
} // namespace BamOtherConstant

#endif // UTILITIES_BAM_BAM_CONSTANT_H_
