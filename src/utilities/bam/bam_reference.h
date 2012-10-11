#ifndef UTILITIES_BAM_BAM_HEADER_H
#define UTILITIES_BAM_BAM_HEADER_H

#include "utilities/bam/SR_BamHeader.h"

class BamReference {
 public:
  BamReference();
  void Init(const SR_BamHeader& bam_header);
  inline const int GetCount() const;
  inline const char* GetName(const int& id) const;

  int             count_no_special;
 private:
  int             count;
  const char**    names;
  const uint32_t* lengths;
  bool isInited;

  BamReference (const BamReference&);
  BamReference& operator=(const BamReference&);
};

inline const int BamReference::GetCount() const { 
  return count;
}

inline const char* BamReference::GetName(const int& id) const {
  if ((id > count) || !isInited)
    return NULL;

  return names[id];
}

#endif // UTILITIES_BAM_BAM_HEADER_H
