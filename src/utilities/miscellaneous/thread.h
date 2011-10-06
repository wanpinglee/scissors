#ifndef UTILITIES_MISCELLANEOUS_THREADDATA_H
#define UTILITIES_MISCELLANEOUS_THREADDATA_H

#include "utilities/bam/SR_BamInStream.h"
#include "utilities/common/SR_Types.h"

struct ThreadData{
  int   id;
  float allowed_clip;
  SR_BamInStream* bam_reader;
  SR_Status bam_status;
};

void StartThreadOrDie (const int& thread_count,
                       SR_BamInStream* bam_reader);

#endif // UTILITIES_MISCELLANEOUS_THREADDATA_H
