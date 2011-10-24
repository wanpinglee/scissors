#include "thread.h"

#include <pthread.h>
#include <stdio.h>

#include <iostream>
#include <vector>

extern "C" {
#include "dataStructures/SR_QueryRegion.h"
#include "utilities/bam/SR_BamPairAux.h"
}

#include "utilities/miscellaneous/aligner.h"

using std::vector;
using std::cout;
using std::endl;

pthread_mutex_t bam_in_mutex;
pthread_mutex_t bam_out_mutex;

// Given a SR_BamListIter containing alignments,
//  reports the chromosome id of the alignments.
// Note: Alignments in the list should be located in the same chromosome
inline void GetChromosomeId(const SR_BamListIter& alignment_list,
                     int* chromosome_id) {
  *chromosome_id = alignment_list->alignment.core.tid;
}

void StoreAlignmentInBam(const vector<bam1_t>& alignments,
                         bamFile* bam_writer) {
  pthread_mutex_lock(&bam_out_mutex);
  for (unsigned int i = 0; i < alignments.size(); ++i) {
    /*
    cout << alignments[i].query_name << "\t" 
         << thread_id << "\t"
         << alignments[i].reference_index << "\t"
	 << alignments[i].reference_begin << "\t"
	 << alignments[i].reference_end << endl;
    */
  }
  pthread_mutex_unlock(&bam_out_mutex);
}

void* RunThread (void* thread_data_) {
  ThreadData *td = (ThreadData*) thread_data_;
  Aligner aligner(td->reference, td->hash_table);

  while (true) { // until bam != SR_OK
    SR_Status bam_status;

    // try to get alignments
    pthread_mutex_lock(&bam_in_mutex);
    bam_status = *(td->bam_status);
    bool terminate = false;
    if (td->alignment_list == NULL) {
      pthread_mutex_lock(&bam_out_mutex);
      cout << "**" << td->id << endl;
      pthread_mutex_unlock(&bam_out_mutex);

      if (bam_status == SR_OK) {
        bam_status = SR_LoadUniquOrphanPairs(td->bam_reader, 
                                             td->id, 
                                             td->allowed_clip);
        *(td->bam_status) = bam_status;
        td->alignment_list = SR_BamInStreamGetIter(td->bam_reader,
                                                   td->id);
      } else {
        terminate = true; // break the while loop
      }
    }

    pthread_mutex_unlock(&bam_in_mutex);
    if (terminate) break;

    if (td->alignment_list != NULL) {
      td->alignments.clear();
      aligner.AlignCandidate(&td->alignment_list, &td->alignments);
      StoreAlignmentInBam(td->alignments, td->bam_writer);
      
      pthread_mutex_lock(&bam_in_mutex);
      SR_BamInStreamClearBuff(td->bam_reader, td->id);
      pthread_mutex_unlock(&bam_in_mutex);
    } else {
      break; // break the while loop
    }

  } // end while

  pthread_exit(NULL);
}

Thread::Thread(const BamReference* bam_reference,
	       const float& allowed_clip,
	       const int& thread_count,
	       FILE* ref_reader,
	       FILE* hash_reader,
	       SR_BamInStream* bam_reader,
	       bamFile*        bam_writer)
    : bam_reference_(bam_reference)
    , allowed_clip_(allowed_clip)
    , thread_count_(thread_count)
    , ref_reader_(ref_reader)
    , hash_reader_(hash_reader)
    , bam_reader_(bam_reader)
    , bam_writer_(bam_writer)
{
  bam_status_ = SR_OK;
  InitThreadData();
  Init();

}

void Thread::Init() {
  reference_        = SR_ReferenceAlloc();
  reference_header_ = SR_RefHeaderAlloc();

  int64_t reference_seal = SR_RefHeaderRead(reference_header_, ref_reader_);
  unsigned char hash_size = 0;
  int64_t hash_seal = SR_InHashTableReadStart(&hash_size, hash_reader_);
  if (reference_seal != hash_seal) {
    printf("ERROR: The reference file is not compatible with the hash table file.\n");
    exit(1);
  }

  hash_table_ = SR_InHashTableAlloc(hash_size);

}

Thread::~Thread() {
  SR_ReferenceFree(reference_);
  SR_InHashTableFree(hash_table_);
  SR_RefHeaderFree(reference_header_);
}

void Thread::InitThreadData() {
  thread_data_.resize(thread_count_);
  for (int i = 0; i < thread_count_; ++i) {
    thread_data_[i].id             = i;
    thread_data_[i].allowed_clip   = allowed_clip_;
    thread_data_[i].bam_reader     = bam_reader_;
    thread_data_[i].alignment_list = NULL;
    thread_data_[i].bam_status     = &bam_status_;
    thread_data_[i].reference      = reference_;
    thread_data_[i].hash_table     = hash_table_;
    thread_data_[i].bam_writer     = bam_writer_;
    thread_data_[i].alignments.clear();
    SR_BamInStreamClearBuff(bam_reader_, i);
  }
}

// Note: we need to load an alignment to know which reference and 
//   hash table that we are gonna load
bool Thread::LoadReference() {
  int thread_id = 0;
  
  while ((thread_data_[0].alignment_list == NULL) &&
        (bam_status_ != SR_EOF)) {
    bam_status_ = SR_LoadUniquOrphanPairs(bam_reader_,
                                          thread_id,
  				          allowed_clip_);
    cout << bam_status_ << endl;
    if (bam_status_ == SR_ERR) { // cannot load alignments from bam
      cout << "ERROR: Cannot load alignments from the input bam." << endl;
      return false;
    }

    thread_data_[0].alignment_list = SR_BamInStreamGetIter(bam_reader_, 
                                                           thread_id);
    if (thread_data_[0].alignment_list == NULL)
      cout << "NULL" << endl;
    else
      cout << "Not NULL" << endl;
  }

  if (thread_data_[0].alignment_list == NULL)
  // this also means bam_status_ != SR_EOF
    return true;
 
  int chromosome_id;
  GetChromosomeId(thread_data_[0].alignment_list, &chromosome_id);
  cout << chromosome_id << endl;

  if (chromosome_id > bam_reference_->GetCount()) {
  // the obtained chr id is invalid
    cout << "ERROR: Reference id is larger than total references." << endl;
    return false;
  }
  
  const char* ref_name = bam_reference_->GetName(chromosome_id);
  int32_t ref_id = SR_RefHeaderGetRefID(reference_header_, ref_name);

  if (ref_id < 0) {  // Cannot find the reference
    printf("ERROR: The reference, %s, is not found in reference file.\n", ref_name);
	exit(1);
  } else {
    SR_ReferenceJump(ref_reader_, reference_header_, ref_id);
    SR_InHashTableJump(hash_reader_, reference_header_, ref_id);
    SR_ReferenceRead(reference_, ref_reader_);
    SR_InHashTableRead(hash_table_, hash_reader_);
  } // end if-else
  return true;
}

bool Thread::Start() {
  vector<pthread_t> threads;
  threads.resize(thread_count_);

  pthread_attr_t attr;

  while (true) { // break when bam_status != SR_OUT_OF_RANGE
    InitThreadData();
    if (!LoadReference()) // load the next first chunk of alignments
      return false;       // and based on the chr id in alignments
		          // load reference and hash table
    
    // register mutex
    pthread_mutex_init(&bam_in_mutex, NULL);
    pthread_mutex_init(&bam_out_mutex, NULL);
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // run threads
    for (int i = 0; i < thread_count_; ++i) {
      int rc = pthread_create(&threads[i], &attr, RunThread, (void*)&thread_data_[i]);
      if (rc) {
        printf("ERROR: Return code from pthread_create is %d.", rc);
        return false;
      } // end if
    } // end for

    // join threads
    for (int i = 0; i < thread_count_; ++i) {
      void* status;
      int rc = pthread_join(threads[i], &status);
      if (rc) {
        printf("ERROR: Return code from pthread_join is %d.", rc);
        return false;
      }
    } // end for

    if (bam_status_ != SR_OUT_OF_RANGE)
      break; // break the while loop
  } // end while

  return true;
}
