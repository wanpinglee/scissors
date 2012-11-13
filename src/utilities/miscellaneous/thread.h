#ifndef UTILITIES_MISCELLANEOUS_THREAD_H
#define UTILITIES_MISCELLANEOUS_THREAD_H

#include <vector>

extern "C" {
#include "outsources/samtools/bam.h"
#include "utilities/bam/bam_reference.h"
#include "utilities/bam/SR_BamInStream.h"
#include "utilities/common/SR_Types.h"
#include "utilities/hashTable/SR_InHashTable.h"
#include "utilities/hashTable/SR_Reference.h"
}

#include "dataStructures/alignment.h"
#include "dataStructures/target_event.h"
#include "dataStructures/target_region.h"
#include "dataStructures/technology.h"
#include "utilities/hashTable/special_hasher.h"
#include "utilities/hashTable/reference_hasher.h"
#include "utilities/miscellaneous/alignment_filter.h"

using std::vector;

class FastaReference;

namespace Scissors {

struct ThreadData{
  int             id;
  float           allowed_clip;
  //int             fragment_length;
  Technology      technology;
  TargetEvent     target_event;
  int             bam_mq_threshold;
  AlignmentFilter alignment_filter;
  TargetRegion    target_region;
  SR_BamInStream*     bam_reader;
  SR_BamInStreamIter  alignment_list;
  SR_Status*          bam_status;
  SR_Reference*       reference;
  SR_InHashTable*     hash_table;
  SR_Reference*       reference_special;
  SR_InHashTable*     hash_table_special;
  SR_RefHeader*       reference_header;
  bamFile*            bam_writer;
  bamFile*            bam_writer_complete_bam;
  vector<Alignment>   alignments;
  vector<bam1_t*>     alignments_bam;
  vector<bam1_t*>     alignments_anchor;
};

class Thread {
 public: 
  Thread(const BamReference* bam_reference,
         const float&           allowed_clip,
         const int&             thread_count,
//	 const int&             fragment_length,
	 const Technology&      technology,
	 const TargetEvent&     target_event,
	 const int&             bam_mq_threshold,
	 const AlignmentFilter& alignment_filter,
	 const TargetRegion&    target_region,
	 const string           special_fasta,
	 FastaReference*        ref_reader,
	 SR_BamInStream* bam_reader,
	 bamFile*        bam_writer,
	 bamFile*        bam_writer_complete_bam);
  ~Thread();
 bool Start();
 private:
  const BamReference*   bam_reference_;
  const float           allowed_clip_;
  const int             thread_count_;
  //const int             fragment_length_;
  const Technology      technology_;
  const TargetEvent     target_event_;
  const int             bam_mq_threshold_;
  const AlignmentFilter alignment_filter_;
  const TargetRegion    target_region_;
  const string    special_fasta_;
  FastaReference* ref_reader_;
  SR_BamInStream* bam_reader_;
  bamFile*        bam_writer_;
  bamFile*        bam_writer_complete_bam_;
  SR_Status       bam_status_;
  //SR_Reference*   reference_;
  //SR_Reference*   reference_special_;
  //SR_InHashTable* hash_table_;
  //SR_InHashTable* hash_table_special_;
  //SR_RefHeader*   reference_header_;
  vector<ThreadData> thread_data_;
  ReferenceHasher ref_hasher_;
  SpecialHasher   sp_hasher_;
  string reference_bases_;

  void Init();
  void InitThreadData();
  void SetReferenceToThread();
  bool LoadReference();
  Thread (const Thread&);
  Thread& operator=(const Thread&);
};
} //namespace
#endif // UTILITIES_MISCELLANEOUS_THREADDATA_H
