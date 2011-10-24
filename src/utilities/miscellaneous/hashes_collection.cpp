#include "hashes_collection.h"

#include <stdio.h>
#include <algorithm>

bool OperatorLength(const BestRegion* r1, const BestRegion* r2) {
  return (r1->length) < (r2->length);
}

void HashesCollection::Init(const BestRegionArray& array) {
  hash_regions_.resize(array.size);
  for (unsigned int i = 0; i < array.size; ++i) 
    hash_regions_[i] = &array.data[i];
}

void HashesCollection::SortByLength() {
  sort(hash_regions_.begin(), hash_regions_.end(), OperatorLength);
}

const BestRegion* HashesCollection::Get (const unsigned int& index) const {
  if (index >= hash_regions_.size()) return NULL;

  return hash_regions_[index];
}

void HashesCollection::Print(void) const {
  for (vector<BestRegion*>::const_iterator ite = hash_regions_.begin();
      ite != hash_regions_.end(); ++ite) {
    printf("%u %u %u %u %u %u\n", (*ite)->refBegins[0], (*ite)->refBegins[1], (*ite)->refBegins[2], (*ite)->queryBegin, (*ite)->length, (*ite)->numPos);
  }
}
