#ifndef UTILITIES_MISCELLANEOUS_HASH_REGION_COLLECTION_H_ 
#define UTILITIES_MISCELLANEOUS_HASH_REGION_COLLECTION_H_

#include <vector>

#include "utilities/hashTable/SR_HashRegionArray.h"

using std::vector;

class HashesCollection {
 public:
  void Init(const BestRegionArray& array);
  void SortByLength(void);
  void Print(void)const;
  const BestRegion* Get (const unsigned int& index) const;
 private:
  vector<BestRegion*> hash_regions_;
}; // end HashRegionCollection

#endif  // UTILITIES_HASH_REGION_COLLECTION_H_
