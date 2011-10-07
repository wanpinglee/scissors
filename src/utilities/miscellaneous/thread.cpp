#include "thread.h"

#include <pthread.h>
#include <stdio.h>

#include <vector>

extern "C" {
#include "utilities/bam/SR_BamPairAux.h"
}

using std::vector;

pthread_mutex_t bam_in_mutex;
pthread_mutex_t bam_out_mutex;

// Given a SR_BamListIter containing alignments,
//  reports the chromosome id of the alignments.
// Note: Alignments in the list should be located in the same chromosome
bool GetChromosomeId(const SR_BamListIter& alignment_list,
                     int* chromosome_id) {
  if (alignment_list == NULL) {
    return false;
  } else {
    *chromosome_id = alignment_list->alignment.core.tid;
    return true;
  }
}

void* RunThread (void* thread_data_) {
  ThreadData *td = (ThreadData*) thread_data_;

  while (true) { // until bam != SR_OK
    SR_Status bam_status;

    // try to get alignments
    pthread_mutex_lock(&bam_in_mutex);
    if (*(td->bam_status) != SR_OK) {
      // the bam is not okay for getting other alignments
      pthread_mutex_unlock(&bam_in_mutex);
      break; // break the while loop
    } else {
      if (td->alignment_list == NULL) {
        bam_status = SR_LoadUniquOrphanPairs(td->bam_reader, 
                                             td->id, 
                                             td->allowed_clip);
        *(td->bam_status) = bam_status;
      }
      pthread_mutex_unlock(&bam_in_mutex);
    } // end if-else

    // no alignment is loaded
    if (bam_status != SR_OK) {
      break; // break the while loop
    } else {
      
    }
  } // end while

  pthread_exit(NULL);
}

Thread::Thread(const BamReference* bam_reference,
	       const float& allowed_clip,
	       const int& thread_count,
	       FILE* ref_reader,
	       FILE* hash_reader,
	       SR_BamInStream* bam_reader)
    : bam_reference_(bam_reference)
    , allowed_clip_(allowed_clip)
    , thread_count_(thread_count)
    , ref_reader_(ref_reader)
    , hash_reader_(hash_reader)
    , bam_reader_(bam_reader)
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
  }
}

// Note: we need to load an alignment to know which reference and 
//   hash table that we are gonna load
bool Thread::LoadReference() {
  int thread_id = 0;
  bam_status_ = SR_LoadUniquOrphanPairs(bam_reader_,
                                          thread_id,
					  allowed_clip_);
  if (bam_status_ != SR_OK) // cannot load alignments from bam
    return false;
  
  thread_data_[0].alignment_list = SR_BamInStreamGetIter(bam_reader_, 
                                                         thread_id);
  int chromosome_id;
  if (!GetChromosomeId(thread_data_[0].alignment_list, &chromosome_id)) 
  // cannot get chr id
    return false;

  if (chromosome_id > bam_reference_->GetCount()) 
  // the obtained chr id is invalid
    return false;
  
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

  InitThreadData();


  while (true) { // break when bam_status != SR_OUT_OF_RANGE
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
