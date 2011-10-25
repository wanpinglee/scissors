#ifndef DATASTRUCTURES_SEARCH_REGION_TYPE_H_
#define DATASTRUCTURES_SEARCH_REGION_TYPE_H_

#include <vector>

#include "technology.h"

using std::vector;

class SearchRegionType {
 public:
  struct RegionType {
    bool upstream; // the reference coordinate of the target region 
                   //   is greater than the anchor
    bool sequence_inverse;    // the sequence in the target region
                              //   is inverse of the reference
    bool sequence_complement; // the sequence in the target region
                              //  is complement of the reference
  };

  // Constructor
  SearchRegionType();
  explicit SearchRegionType(const Technology& technology);

  // Destructor
  //~SearchRegionType();

  bool GetNextRegionType(const bool is_anchor_forward, RegionType* region_type);
  bool SetCurrentTypeSuccess(const bool is_anchor_forward);
  void ResetRegionTypeList(void);
  inline void RewindRegionTypeList(void);
  inline void SetTechnology(const Technology& technology);

 private:
  Technology technology_;
  vector<RegionType> forward_anchor_type_vector_;
  vector<RegionType> reverse_anchor_type_vector_;
  vector<int> forward_anchor_type_preference_;
  vector<int> reverse_anchor_type_preference_;
  int current_forward_anchor_type_preference_; // id for forward_anchor_type_preference_
  int current_reverse_anchor_type_preference_; // id for reverse_anchor_type_preference_

  vector<int> forward_anchor_type_count_; // corresponds to forward_anchor_type_vector_
  vector<int> reverse_anchor_type_count_; // corresponds to reverse_anchor_type_vector_

  bool has_gotten_forward_type_;
  bool has_gotten_reverse_type_;

  void Init(void);
  SearchRegionType (const SearchRegionType&);
  SearchRegionType& operator= (const SearchRegionType&);

};


// inline functions
inline void SearchRegionType::RewindRegionTypeList(void) {
  current_forward_anchor_type_preference_ = 0;
  current_reverse_anchor_type_preference_ = 0;
  has_gotten_forward_type_ = false;
  has_gotten_reverse_type_ = false;
}

inline void SearchRegionType::SetTechnology(const Technology& technology) {
  technology_ = technology;
  ResetRegionTypeList();
}

#endif // DATASTRUCTURES_SEARCH_REGION_TYPE_H_
