#include "special_hasher.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

extern "C" {
#include "SR_InHashTable.h"
#include "SR_OutHashTable.h"
#include "SR_Reference.h"
}

using std::string;

namespace {
void ConvertHashTableOutToIn(const SR_OutHashTable* out, SR_InHashTable* in) {
  in->id = out->id;

  // indices are already allocated in SR_InHashTableAlloc
  // copy indices from out hash table
  uint32_t index = 1;
  for (unsigned int i = 0; i != out->numHashes; ++i) {
    in->indices[i] = index;
    index += (out->hashPosTable)[i].size;
  }

  // total number of hash positons that are found in the references
  in->numPos = out->numPos;

  free(in->hashPos);
  in->hashPos = (uint32_t*) malloc(sizeof(uint32_t) * in->numPos);
  // sanity checker
  if (in->hashPos == NULL)
    fprintf(stderr, "ERROR: Not enough memory for the storage \
            of hash positions in the hash table object.\n");
  uint32_t* ptr = in->hashPos;
  for (unsigned int i = 0; i != out->numHashes; ++i) {
    uint32_t size = (out->hashPosTable)[i].size;
    memcpy(ptr, (out->hashPosTable)[i].data, sizeof(uint32_t) * size);
    ptr += size;
  }
}
}

SpecialHasher::SpecialHasher()
    : fasta_("")
    , reference_header_(NULL)
    , references_(NULL)
    , hash_table_(NULL)
    , hash_size_(7){
  Init();
}

SpecialHasher::SpecialHasher(const char* fasta, const int& hash_size)
    : fasta_(fasta)
    , reference_header_(NULL)
    , references_(NULL)
    , hash_table_(NULL)
    , hash_size_(hash_size){
  Init();
}

SpecialHasher::~SpecialHasher() {
  // SR_RefHeaderFree also frees memory 
  //   that is allocated by SR_SpecialRefInfoAlloc
  SR_RefHeaderFree(reference_header_);
  SR_ReferenceFree(references_);
  SR_InHashTableFree(hash_table_);
}

void SpecialHasher::Init() {
  // allocate memory for the header and special info
  reference_header_ = SR_RefHeaderAlloc(1,0);
  const int initial_reference_number = 30;
  reference_header_->pSpecialRefInfo = 
    SR_SpecialRefInfoAlloc(initial_reference_number);

  references_ = SR_ReferenceAlloc();
}

bool SpecialHasher::Load() {
  if (fasta_.empty()) {
    fprintf(stderr, "ERROR: Please set fasta filename before loading.\n");
    return false;
  }
  // load the special references
  FILE* input = fopen(fasta_.c_str(), "r");
  if (input == NULL) {
    fprintf(stderr, "ERROR: The file (%s) cannot be opened.\n", fasta_.c_str());
  } else {
    SR_SpecialRefLoad(references_, reference_header_, input);
  }
  fclose(input);

  // index every possible hash position in the current chromosome
  // and write the results into hash position index file and hash position file
  SR_OutHashTable* out_hash_table = SR_OutHashTableAlloc(hash_size_);
  SR_OutHashTableLoad(out_hash_table, references_->sequence, 
                      references_->seqLen, references_->id);

  hash_table_ = SR_InHashTableAlloc(hash_size_);
  ConvertHashTableOutToIn(out_hash_table, hash_table_);

  SR_OutHashTableFree(out_hash_table);

  return true;
}
