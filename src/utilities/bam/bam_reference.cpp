#include "bam_reference.h"

#include "utilities/bam/SR_BamHeader.h"

BamReference::BamReference()
    : count(0)
    , names(NULL)
    , lengths(NULL)
    , isInited(false)
{}

void BamReference::Init(const SR_BamHeader& bam_header) {
  count   = SR_BamHeaderGetRefNum(&bam_header);
  names   = SR_BamHeaderGetRefNames(&bam_header);
  lengths = SR_BamHeaderGetRefLens(&bam_header);

  isInited = true;
}
