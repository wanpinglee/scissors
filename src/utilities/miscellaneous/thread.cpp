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

void StoreAlignmentInBam(const vector<bam1_t*>& alignments_bam,
                         bamFile* bam_writer) {
  pthread_mutex_lock(&bam_out_mutex);
  for (unsigned int i = 0; i < alignments_bam.size(); ++i) {
    bam_write1(*bam_writer, alignments_bam[i]);
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

void FreeAlignmentBam(vector<bam1_t*>* als_bam) {
  
  for (unsigned int i = 0; i < als_bam->size(); ++i) {
    bam1_t* ptr = (*als_bam)[i];
    bam_destroy1(ptr);
  }
  

  als_bam->clear();
}

/*
void ConvertAlignments(const vector<Alignment>& als,
                       vector<bam1_t>* als_bam) {
  als_bam.resize(als.size());

  for (unsigned int i = 0; i < als.size(); ++i) {
    ConvertAlignmentToBam1(als[i], &als_bam[i]);
  }

}
*/

void* RunThread (void* thread_data_) {
  ThreadData *td = (ThreadData*) thread_data_;
  Aligner aligner(td->reference, td->hash_table, td->reference_special, 
                  td->hash_table_special, td->reference_header, td->fragment_length);

  while (true) { // until bam != SR_OK
    SR_Status bam_status;

    // try to get alignments from bam_reader
    pthread_mutex_lock(&bam_in_mutex);
    bam_status = *(td->bam_status);
    bool terminate = false;
    if (td->alignment_list.pBamNode == NULL) {
      if (bam_status == SR_OK) {
        // TODO @WP: make sure each field of Jiantao
        bam_status = SR_LoadAlgnPairs(td->bam_reader,
					     NULL, 
					     // the pointer to frag length 
					     // distribution; NULL means
					     // we don't want to load it
                                             td->id, 
                                             td->allowed_clip,
					     0.1, //maxMismatchRate
					     1 // min mapping quality
					     );
        *(td->bam_status) = bam_status;
        SR_BamInStreamSetIter(&td->alignment_list,
	                      td->bam_reader,
                              td->id);
      } else {
        terminate = true; // break the while loop
      }
    }

    pthread_mutex_unlock(&bam_in_mutex);
    if (terminate) break;

    if (td->alignment_list.pBamNode != NULL) {
      td->alignments.clear();
      aligner.AlignCandidate(td->detect_special, &td->alignment_list, &td->alignments_bam);
      StoreAlignmentInBam(td->alignments_bam, td->bam_writer);
      FreeAlignmentBam(&td->alignments_bam);
      
      pthread_mutex_lock(&bam_in_mutex);
      SR_BamInStreamClearRetList(td->bam_reader, td->id);
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
	       const int& fragment_length,
	       const bool& detect_special,
	       FILE* ref_reader,
	       FILE* hash_reader,
	       SR_BamInStream* bam_reader,
	       bamFile*        bam_writer)
    : bam_reference_(bam_reference)
    , allowed_clip_(allowed_clip)
    , thread_count_(thread_count)
    , fragment_length_(fragment_length)
    , detect_special_(detect_special)
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
  int64_t reference_seal;
  // load the reference header
  reference_header_ = SR_RefHeaderRead(&reference_seal, ref_reader_);
  
  // load the header in hash table file
  unsigned char hash_size = 0;
  int64_t hash_seal = SR_InHashTableReadStart(&hash_size, hash_reader_);
  
  // check the compatibility between reference file and the hash table file
  if (reference_seal != hash_seal) {
    printf("ERROR: The reference file is not compatible with the hash table file.\n");
    exit(1);
  }

  // load special references and their hash tables if necessary
  if (detect_special_) {
    // load special references
    reference_special_ = SR_ReferenceAlloc();
    if (SR_SpecialRefRead(reference_special_, reference_header_, ref_reader_) == SR_ERR)
      printf("ERROR: Cannot read special sequence.\n");
    // load special hash tables
    hash_table_special_ = SR_InHashTableAlloc(hash_size);
    if (SR_InHashTableReadSpecial(hash_table_special_, reference_header_, hash_reader_) == SR_ERR)
      printf("ERROR: Cannot read special hash table.\n");
  }

  hash_table_ = SR_InHashTableAlloc(hash_size);
  reference_  = SR_ReferenceAlloc();

}

Thread::~Thread() {
  SR_ReferenceFree(reference_);
  SR_InHashTableFree(hash_table_);
  SR_RefHeaderFree(reference_header_);
  if (detect_special_) {
    SR_ReferenceFree(reference_special_);
    SR_InHashTableFree(hash_table_special_);
  }
}

void Thread::InitThreadData() {
  thread_data_.resize(thread_count_);
  for (int i = 0; i < thread_count_; ++i) {
    thread_data_[i].id                       = i;
    thread_data_[i].allowed_clip             = allowed_clip_;
    thread_data_[i].fragment_length          = fragment_length_;
    thread_data_[i].detect_special           = detect_special_;
    thread_data_[i].bam_reader               = bam_reader_;
    thread_data_[i].alignment_list.pBamNode  = NULL;
    thread_data_[i].alignment_list.pAlgnType = NULL;
    thread_data_[i].bam_status               = &bam_status_;
    thread_data_[i].reference                = reference_;
    thread_data_[i].hash_table               = hash_table_;
    thread_data_[i].reference_special        = reference_special_;
    thread_data_[i].hash_table_special       = hash_table_special_;
    thread_data_[i].reference_header         = reference_header_;
    thread_data_[i].bam_writer               = bam_writer_;
    thread_data_[i].alignments.clear();
    FreeAlignmentBam(&thread_data_[i].alignments_bam);
    SR_BamInStreamClearRetList(bam_reader_, i);
  }
}

// Note: we need to load an alignment to know which reference and 
//   hash table that we are gonna load
bool Thread::LoadReference() {
  int thread_id = 0;
  
  while ((thread_data_[0].alignment_list.pBamNode == NULL) &&
        (bam_status_ != SR_EOF)) {
    bam_status_ = SR_LoadAlgnPairs(bam_reader_,
                                      NULL,
				      // the pointer to frag length
				      // distribution; NULL means
				      // we don't want to load it
                                      thread_id,
  				      allowed_clip_,
				      0.1, // maxMismatchRate
				      1); // min mapping quality
    if (bam_status_ == SR_ERR) { // cannot load alignments from bam
      cout << "ERROR: Cannot load alignments from the input bam." << endl;
      return false;
    }

    //cout << bam_status_ << endl;
    if (bam_status_ == SR_OUT_OF_RANGE)
      continue;
    else
    SR_BamInStreamSetIter(&thread_data_[0].alignment_list,
                          bam_reader_, 
                          thread_id);
  }

  if (thread_data_[0].alignment_list.pBamNode == NULL) {
  // this also means bam_status_ != SR_EOF
    return true;
  }
 
  int chromosome_id;
  GetChromosomeId(thread_data_[0].alignment_list.pBamNode, &chromosome_id);

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
