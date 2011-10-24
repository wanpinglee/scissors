#ifndef UTILITIES_MISCELLANEOUS_THREAD_H
#define UTILITIES_MISCELLANEOUS_THREAD_H

#include <vector>

extern "C" {
#include "utilities/bam/bam_reference.h"
#include "utilities/bam/SR_BamInStream.h"
#include "utilities/common/SR_Types.h"
#include "utilities/hashTable/SR_InHashTable.h"
#include "utilities/hashTable/SR_Reference.h"
}

#include "utilities/bam/bam_alignment.h"
#include "utilities/bam/bam_writer.h"

using std::vector;

struct ThreadData{
  int   id;
  float allowed_clip;
  SR_BamInStream* bam_reader;
  SR_BamListIter  alignment_list;
  SR_Status*      bam_status;
  SR_Reference*   reference;
  SR_InHashTable* hash_table;
  BamWriter*      bam_writer;
  vector<BamAlignment> alignments;
};

class Thread {
 public: 
  Thread(const BamReference* bam_reference,
         const float&    allowed_clip,
         const int&      thread_count,
	 FILE*           ref_reader,
         FILE*           hash_reader,
	 SR_BamInStream* bam_reader,
	 BamWriter*      bam_writer);
  ~Thread();
 bool Start();
 private:
  const BamReference*   bam_reference_;
  const float           allowed_clip_;
  const int             thread_count_;
  FILE*           ref_reader_;
  FILE*           hash_reader_;
  SR_BamInStream* bam_reader_;
  BamWriter*      bam_writer_;
  SR_Status       bam_status_;
  SR_Reference*   reference_;
  SR_InHashTable* hash_table_;
  SR_RefHeader*   reference_header_;
  vector<ThreadData> thread_data_;

  void Init();
  void InitThreadData();
  bool LoadReference();
  Thread (const Thread&);
  Thread& operator=(const Thread&);
};

#endif // UTILITIES_MISCELLANEOUS_THREADDATA_H
