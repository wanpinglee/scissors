#ifndef DATASTRUCTURES_SEARCH_REGION_TYPE_H_
#define DATASTRUCTURES_SEARCH_REGION_TYPE_H_

#include <vector>

#include "technology.h"

using std::vector;

namespace Scissors {
class SearchRegionType {
 public:
  struct RegionType {
    bool upstream; // The reference coordinate of the target region 
                   // is greater than the anchor.
    bool sequence_inverse;    // The sequence in the target region
                              // is inverse of the reference.
    bool sequence_complement; // The sequence in the target region
                              // is complement of the reference.
  };

  // Constructor
  SearchRegionType();
  SearchRegionType(const Technology& technology, const bool& mate1);

  // Destructor
  //~SearchRegionType();

  bool GetNextRegionType(const bool is_anchor_forward, RegionType* region_type);
  bool GetStandardType(const bool is_anchor_forward, RegionType* region_type);
  bool SetCurrentTypeSuccess(const bool is_anchor_forward);
  void ResetRegionTypeList(void);
  inline void RewindRegionTypeList(void);
  inline void SetTechnologyAndAnchorMate1(const Technology& technology, const bool& mate1);
  inline void SetTechnology(const Technology& technology);
  inline void SetAnchorMate1(const bool& mate1);

 private:
  Technology technology_;
  bool mate1_;
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

  static const int types_ = 6;

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

inline void SearchRegionType::SetTechnologyAndAnchorMate1(const Technology& technology, const bool& mate1) {
  technology_ = technology;
  mate1_ = mate1;
  Init();
}

inline void SearchRegionType::SetTechnology(const Technology& technology) {
  technology_ = technology;
  Init();
}

inline void SearchRegionType::SetAnchorMate1(const bool& mate1) {
  mate1_ = mate1;
  Init();
}
} // namespace
#endif // DATASTRUCTURES_SEARCH_REGION_TYPE_H_
