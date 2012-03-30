#include <stdlib.h>
#include <string.h>

extern "C" {
#include "SR_InHashTable.h"
#include "SR_OutHashTable.h"
}

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
