#ifndef _ALIGNER_H_
#define _ALIGNER_H_

extern "C" {
#include "dataStructures/SR_QueryRegion.h"
#include "outsources/samtools/bam.h"
#include "utilities/bam/SR_BamInStream.h"
#include "utilities/hashTable/SR_HashRegionTable.h"
#include "utilities/hashTable/SR_InHashTable.h"
#include "utilities/hashTable/SR_Reference.h"
}

#include "dataStructures/alignment.h"
#include "dataStructures/anchor_region.h"
#include "dataStructures/search_region_type.h"

class Aligner {
 public:
  Aligner(const SR_Reference* reference, 
          const SR_InHashTable* hash_table);
  ~Aligner();
  void AlignCandidate(SR_BamListIter* al_ite, vector<bam1_t*>* alignments);
 private:
  SearchRegionType search_region_type_;
  AnchorRegion     anchor_region_;

  const SR_Reference*   reference_;
  const SR_InHashTable* hash_table_;
  SR_QueryRegion*       query_region_;
  HashRegionTable*      hashes_;
  SR_SearchArgs         hash_length_;

  void LoadRegionType(const bam1_t& anchor);
  inline const char* GetSequence(const size_t& start) const;

  Aligner (const Aligner&);
  Aligner& operator= (const Aligner&);
};


inline const char* Aligner::GetSequence(const size_t& start) const {
  if (start >= reference_->seqLen)
    return NULL;
  else
    return (reference_->sequence + start);
}

#endif // _ALIGNER_H_
