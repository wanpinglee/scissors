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
#include "utilities/miscellaneous/alignment_filter.h"

using std::vector;

struct ThreadData{
  int   id;
  float allowed_clip;
  int   fragment_length;
  bool  detect_special;
  AlignmentFilter alignment_filter;
  SR_BamInStream*     bam_reader;
  SR_BamInStreamIter  alignment_list;
  SR_Status*          bam_status;
  SR_Reference*       reference;
  SR_InHashTable*     hash_table;
  SR_Reference*       reference_special;
  SR_InHashTable*     hash_table_special;
  SR_RefHeader*       reference_header;
  bamFile*            bam_writer;
  vector<Alignment>   alignments;
  vector<bam1_t*>     alignments_bam;
};

class Thread {
 public: 
  Thread(const BamReference* bam_reference,
         const float&    allowed_clip,
         const int&      thread_count,
	 const int&      fragment_length,
	 const bool&     detect_special,
	 const AlignmentFilter& alignment_filter,
	 FILE*           ref_reader,
         FILE*           hash_reader,
	 SR_BamInStream* bam_reader,
	 bamFile*        bam_writer);
  ~Thread();
 bool Start();
 private:
  const BamReference*   bam_reference_;
  const float           allowed_clip_;
  const int             thread_count_;
  const int             fragment_length_;
  const bool            detect_special_;
  const AlignmentFilter alignment_filter_;
  FILE*           ref_reader_;
  FILE*           hash_reader_;
  SR_BamInStream* bam_reader_;
  bamFile*        bam_writer_;
  SR_Status       bam_status_;
  SR_Reference*   reference_;
  SR_Reference*   reference_special_;
  SR_InHashTable* hash_table_;
  SR_InHashTable* hash_table_special_;
  SR_RefHeader*   reference_header_;
  vector<ThreadData> thread_data_;

  void Init();
  void InitThreadData();
  bool LoadReference();
  Thread (const Thread&);
  Thread& operator=(const Thread&);
};

#endif // UTILITIES_MISCELLANEOUS_THREADDATA_H
