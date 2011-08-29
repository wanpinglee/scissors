#ifndef DATASTRUCTURES_SEARCH_REGION_TYPE_H_
#define DATASTRUCTURES_SEARCH_REGION_TYPE_H_

#include <list>

#include "technology.h"

using std::list;

struct RegionType {
  bool upstream;
  bool sequence_reverse;
  bool sequence_complement;
};

class SearchRegionType {
 public:
  // Constructor
  SearchRegionType();
  explicit SearchRegionType( const Technology technology );

  // Destructor
  ~SearchRegionType();

  void ResetRegionTypeList();
  void RewindRegionTypeList();
  bool GetNextRegionType(const bool is_anchor_forward, RegionType* region_type);

 private:
  const Technology technology_;
  list<RegionType> forward_anchor_region_type_list_;
  list<RegionType> reverse_anchor_region_type_list_;
  list<RegionType>::iterator forward_anchor_region_type_list_ptr_;
  list<RegionType>::iterator reverse_anchor_region_type_list_ptr_;

};

#endif // DATASTRUCTURES_SEARCH_REGION_TYPE_H_
