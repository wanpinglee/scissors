#include "search_region_type.h" 

// for illumina forward anchors
const static RegionType kRegionType1 = {true, true, true};
const static RegionType kRegionType2 = {false, true, true};
const static RegionType kRegionType3 = {true, false, true};
const static RegionType kRegionType4 = {false, false, true};

// for illumina reversed complement anchors
const static RegionType kRegionType5 = {false, false, false};
const static RegionType kRegionType6 = {true, false, false};
const static RegionType kRegionType7 = {false, true, false};
const static RegionType kRegionType8 = {true, true, false};


SearchRegionType::SearchRegionType()
    : technology_(ILLUMINA) // the default tech is ILLUMINA
{
  forward_anchor_region_type_list_.resize(4);
  reverse_anchor_region_type_list_.resize(4);
  ResetRegionTypeList();
}

SearchRegionType::SearchRegionType(const Technology technology)
    : technology_(technology)
{
  forward_anchor_region_type_list_.resize(4);
  reverse_anchor_region_type_list_.resize(4);
  ResetRegionTypeList();
}

bool SearchRegionType::LoadNextRegionType(const bool is_anchor_forward, 
    RegionType* region_type) {

  list<RegionType>::iterator ptr = is_anchor_forward ? forward_anchor_region_type_list_ptr_ 
      : reverse_anchor_region_type_list_ptr_;

  if (is_anchor_forward) {
    if (forward_anchor_region_type_list_ptr_ == forward_anchor_region_type_list_.end()) return false;
    *region_type = *forward_anchor_region_type_list_ptr_;
    ++forward_anchor_region_type_list_ptr_;
  } else {
    if (reverse_anchor_region_type_list_ptr_ == reverse_anchor_region_type_list_.end()) return false;
    *region_type = *reverse_anchor_region_type_list_ptr_;
    ++reverse_anchor_region_type_list_ptr_;
  }

  return true;
}

void SearchRegionType::RewindRegionTypeList(){
  forward_anchor_region_type_list_ptr_ = forward_anchor_region_type_list_.begin();
  reverse_anchor_region_type_list_ptr_ = reverse_anchor_region_type_list_.begin();
}

void SearchRegionType::ResetRegionTypeList() {
  forward_anchor_region_type_list_ptr_ = forward_anchor_region_type_list_.begin();
  reverse_anchor_region_type_list_ptr_ = reverse_anchor_region_type_list_.begin();
  
  switch (technology_) {
    case ILLUMINA: {
      // for forward anchor region type list
      *forward_anchor_region_type_list_ptr_ = kRegionType1; 
      ++forward_anchor_region_type_list_ptr_;
      *forward_anchor_region_type_list_ptr_ = kRegionType2; 
      ++forward_anchor_region_type_list_ptr_;
      *forward_anchor_region_type_list_ptr_ = kRegionType3; 
      ++forward_anchor_region_type_list_ptr_;
      *forward_anchor_region_type_list_ptr_ = kRegionType4;
      // for reverse anchor region type list
      *reverse_anchor_region_type_list_ptr_ = kRegionType5;
      ++reverse_anchor_region_type_list_ptr_;
      *reverse_anchor_region_type_list_ptr_ = kRegionType6;
      ++reverse_anchor_region_type_list_ptr_;
      *reverse_anchor_region_type_list_ptr_ = kRegionType7;
      ++reverse_anchor_region_type_list_ptr_;
      *reverse_anchor_region_type_list_ptr_ = kRegionType8;
      break;
    } // ILLUMINA
    // TODO(WP): assign the types for LS454 reads
    case LS454: {
      break;
    }

    default: {
    }
  }
}
