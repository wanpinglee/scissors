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

void* RunThread (void* thread_data) {
  ThreadData *td = (ThreadData*) thread_data;
  int thread_id = 0;

  // get thread id
  pthread_mutex_lock(&bam_in_mutex);
  thread_id = td->id;
  ++(td->id);
  pthread_mutex_unlock(&bam_in_mutex);

  while (true) { // until bam != SR_OK
    SR_Status bam_status;

    // try to get alignments
    pthread_mutex_lock(&bam_in_mutex);
    if (td->bam_status != SR_OK) {
      // the bam is not okay for getting other alignments
      pthread_mutex_unlock(&bam_in_mutex);
      break; // break the while loop
    } else {
      bam_status = 
        SR_LoadUniquOrphanPairs(td->bam_reader, thread_id, td->allowed_clip);
      td->bam_status = bam_status;
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

void StartThreadOrDie (const int& thread_count,
                       SR_BamInStream* bam_reader) {
  vector<pthread_t> threads;
  threads.resize(thread_count);

  // prepare thread attr
  pthread_attr_t  attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  // register mutex
  pthread_mutex_init(&bam_in_mutex, NULL);
  pthread_mutex_init(&bam_out_mutex, NULL);

  ThreadData td;
  td.id = 0;
  td.bam_reader = bam_reader;

  while (true) { // break when bam_status != SR_OUT_OF_RANGE
  // run threads
  for (int i = 0; i < thread_count; ++i) {
    int rc = pthread_create(&threads[i], &attr, RunThread, (void*)&td);
    if (rc) {
      printf("ERROR: Return code from pthread_create is %d.", rc);
      exit(1);
    } // end if
  } // end for

  // join threads
  for (int i = 0; i < thread_count; ++i) {
    void* status;
    int rc = pthread_join(threads[i], &status);
    if (rc) {
      printf("ERROR: Return code from pthread_join is %d.", rc);
      exit(1);
    }
  } // end for

  if (td.bam_status != SR_OUT_OF_RANGE)
    break; // break the while loop
  } // end while
}
