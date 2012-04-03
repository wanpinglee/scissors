#include "reference_hasher.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <string>

extern "C" {
#include "SR_InHashTable.h"
#include "SR_OutHashTable.h"
#include "SR_Reference.h"
}

#include "ConvertHashTableOutToIn.h"

using std::string;

ReferenceHasher::ReferenceHasher(void)
    : references_(NULL)
    , hash_table_(NULL)
    , hash_size_(7)
    , is_loaded_(false){
  Init();
}

ReferenceHasher::ReferenceHasher(const char* sequence)
    : references_(NULL)
    , hash_table_(NULL)
    , hash_size_(7)
    , is_loaded_(false){
  Init();
  SetSequence(sequence);
}

ReferenceHasher::~ReferenceHasher(void) {
  delete references_;
  SR_InHashTableFree(hash_table_);
}

void ReferenceHasher::Init(void) {
  references_ = new SR_Reference;
  references_->sequence = NULL;
  references_->id       = 0;
  references_->seqLen   = 0;
  references_->seqCap   = 0;
}

void ReferenceHasher::SetSequence(const char* sequence) {
  references_->sequence = (char*)sequence;
  references_->seqLen   = strlen(sequence);
}

bool ReferenceHasher::Load(void) {
  if ((references_->sequence == NULL) || (references_->seqLen == 0)) {
    fprintf(stderr, "ERROR: Please set the reference sequence before loading.\n");
    fprintf(stderr, "       The reference length is %u.\n", references_->seqLen);
    return false;
  }

  // index every possible hash position in the current chromosome
  // and write the results into hash position index file and hash position file
  SR_OutHashTable* out_hash_table = SR_OutHashTableAlloc(hash_size_);
  SR_OutHashTableLoad(out_hash_table, references_->sequence, 
                      references_->seqLen, references_->id);

  hash_table_ = SR_InHashTableAlloc(hash_size_);
  ConvertHashTableOutToIn(out_hash_table, hash_table_);

  SR_OutHashTableFree(out_hash_table);

  is_loaded_ = true;
  return true;
}
