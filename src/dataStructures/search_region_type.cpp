#include "search_region_type.h" 

#include <stdio.h>

#include <algorithm>

namespace Scissors {
// struct RegionType {bool upstream; bool sequence_inverse; bool sequence_complement;};
const static SearchRegionType::RegionType kRegionType1 = {true,  true,  true};
const static SearchRegionType::RegionType kRegionType2 = {false, true,  true};
const static SearchRegionType::RegionType kRegionType3 = {true,  false, true};
const static SearchRegionType::RegionType kRegionType4 = {false, false, true};

const static SearchRegionType::RegionType kRegionType5 = {false, false, false};
const static SearchRegionType::RegionType kRegionType6 = {true,  false, false};
const static SearchRegionType::RegionType kRegionType7 = {false, true,  false};
const static SearchRegionType::RegionType kRegionType8 = {true,  true,  false};

SearchRegionType::SearchRegionType()
    : technology_(TECH_ILLUMINA) // the default tech is ILLUMINA
    , mate1_(true)
{
  forward_anchor_type_vector_.resize(types_);
  reverse_anchor_type_vector_.resize(types_);
  forward_anchor_type_preference_.resize(types_);
  reverse_anchor_type_preference_.resize(types_);
  forward_anchor_type_count_.resize(types_);
  reverse_anchor_type_count_.resize(types_);
  Init();
}

SearchRegionType::SearchRegionType(const Technology& technology, const bool& mate1)
    : technology_(technology)
    , mate1_(mate1)
{
  forward_anchor_type_vector_.resize(types_);
  reverse_anchor_type_vector_.resize(types_);
  forward_anchor_type_preference_.resize(types_);
  reverse_anchor_type_preference_.resize(types_);
  forward_anchor_type_count_.resize(types_);
  reverse_anchor_type_count_.resize(types_);
  Init();
}

void SearchRegionType::Init() {
  switch (technology_) {
    case TECH_ILLUMINA: {
      // for forward anchor region type list
      forward_anchor_type_vector_[0] = kRegionType1; // (TTT) normal; pair-order, reverse, and complement
      forward_anchor_type_vector_[1] = kRegionType2; // (FTT) pair-disorder, reverse, and complement
      forward_anchor_type_vector_[2] = kRegionType6; // (TFF) pair-order, not reverse, and not complement
      forward_anchor_type_vector_[3] = kRegionType5; // (FFF) pair-disorder, not reverse, and not complement
      forward_anchor_type_vector_[4] = kRegionType3; // (TFT) pair-order, not reverse, and complement
      forward_anchor_type_vector_[5] = kRegionType4; // (FFT) pair-disorder, not reverse, and complement
      // for reverse anchor region type list
      reverse_anchor_type_vector_[0] = kRegionType5; // (FFF) normal; pair-order, not reverse, and not complement
      reverse_anchor_type_vector_[1] = kRegionType6; // (TFF) pair-disorder, not reverse, and not complement
      reverse_anchor_type_vector_[2] = kRegionType2; // (FTT) pair-order, reverse, and complement
      reverse_anchor_type_vector_[3] = kRegionType1; // (TTT) pair-disorder, reverse, and complement
      reverse_anchor_type_vector_[4] = kRegionType7; // (FTF) pair-order, reverse, and not complement
      reverse_anchor_type_vector_[5] = kRegionType8; // (TTF) pair-disorder, reverse, and not complement
      break;
    } // ILLUMINA
    
    case TECH_SOLID: {
      if (mate1_) {
        // for forward anchor region type list
        // ----->(m1)   ----->(m2)
        forward_anchor_type_vector_[0] = kRegionType6; // (TFF)
        forward_anchor_type_vector_[1] = kRegionType5; // (FFF)
        forward_anchor_type_vector_[2] = kRegionType1; // (TTT)
        forward_anchor_type_vector_[3] = kRegionType2; // (FTT)
        forward_anchor_type_vector_[4] = kRegionType8; // (TTF)
        forward_anchor_type_vector_[5] = kRegionType7; // (FTF)
        // for reverse anchor region type list
        // <-----(m2)   <-----(m1)
        reverse_anchor_type_vector_[0] = kRegionType2; // (FTT)
        reverse_anchor_type_vector_[1] = kRegionType1; // (TTT)
        reverse_anchor_type_vector_[2] = kRegionType5; // (FFF)
        reverse_anchor_type_vector_[3] = kRegionType6; // (TFF)
        reverse_anchor_type_vector_[4] = kRegionType4; // (FFT)
        reverse_anchor_type_vector_[5] = kRegionType3; // (TFT)
      } else {
        // for forward anchor region type list
        // ----->(m1)   ----->(m2)
        forward_anchor_type_vector_[0] = kRegionType5; // (FFF)
        forward_anchor_type_vector_[1] = kRegionType6; // (TFF)
        forward_anchor_type_vector_[2] = kRegionType2; // (FTT)
        forward_anchor_type_vector_[3] = kRegionType1; // (TTT)
        forward_anchor_type_vector_[4] = kRegionType7; // (FTF)
        forward_anchor_type_vector_[5] = kRegionType8; // (TTF)
        // for reverse anchor region type list
        // <-----(m2)   <-----(m1)
        reverse_anchor_type_vector_[0] = kRegionType1; // (TTT)
        reverse_anchor_type_vector_[1] = kRegionType2; // (FTT)
        reverse_anchor_type_vector_[2] = kRegionType6; // (TFF)
        reverse_anchor_type_vector_[3] = kRegionType5; // (FFF)
        reverse_anchor_type_vector_[4] = kRegionType3; // (TFT)
	reverse_anchor_type_vector_[5] = kRegionType4; // (FFT)
      }
      break;
    }
    
    // TODO(WP): assign the types for LS454 reads
    case TECH_454: {
      if (mate1_) {
        // for forward anchor region type list
	// ----->(m2)   ----->(m1)
        forward_anchor_type_vector_[0] = kRegionType5; // (FFF)
        forward_anchor_type_vector_[1] = kRegionType6; // (TFF)
        forward_anchor_type_vector_[2] = kRegionType2; // (FTT)
        forward_anchor_type_vector_[3] = kRegionType1; // (TTT)
        forward_anchor_type_vector_[4] = kRegionType7; // (FTF)
        forward_anchor_type_vector_[5] = kRegionType8; // (TTF)
        // for reverse anchor region type list
	// <-----(m1)   <-----(m2)
        reverse_anchor_type_vector_[0] = kRegionType1; // (TTT)
        reverse_anchor_type_vector_[1] = kRegionType2; // (FTT)
        reverse_anchor_type_vector_[2] = kRegionType6; // (TFF)
        reverse_anchor_type_vector_[3] = kRegionType5; // (FFF)
        reverse_anchor_type_vector_[4] = kRegionType3; // (TFT)
        reverse_anchor_type_vector_[5] = kRegionType4; // (FFT)
      } else {
        // for forward anchor region type list
	// ----->(m2)   ----->(m1)
        forward_anchor_type_vector_[0] = kRegionType6; // (TFF)
        forward_anchor_type_vector_[1] = kRegionType5; // (FFF)
        forward_anchor_type_vector_[2] = kRegionType1; // (TTT)
        forward_anchor_type_vector_[3] = kRegionType2; // (FTT)
        forward_anchor_type_vector_[4] = kRegionType8; // (TTF)
        forward_anchor_type_vector_[5] = kRegionType7; // (FTF)
        // for reverse anchor region type list
	// <-----(m1)   <-----(m2)
        reverse_anchor_type_vector_[0] = kRegionType2; // (FTT)
        reverse_anchor_type_vector_[1] = kRegionType1; // (TTT)
        reverse_anchor_type_vector_[2] = kRegionType5; // (FFF)
        reverse_anchor_type_vector_[3] = kRegionType6; // (TFF)
        reverse_anchor_type_vector_[4] = kRegionType4; // (FFT)
        reverse_anchor_type_vector_[5] = kRegionType3; // (TFT)
      }
      break;
    } // LS454

    default: {
    } // default
  } // switch


  ResetRegionTypeList();
  RewindRegionTypeList();
}


// Use this function when the current type successfully supports 
// split-read alignment.
// The function will rewind the list since the current is successful,
// other types should not be applied.
bool SearchRegionType::SetCurrentTypeSuccess(const bool is_anchor_forward){
  // No type is gotten, so shouldn't set success.
  if (is_anchor_forward && !has_gotten_forward_type_)  return false;
  if (!is_anchor_forward && !has_gotten_reverse_type_) return false;

  // increases the current counter and update the list if necessary
  // if the current counter is greater than the previous
  // then swap the current and the previous one
  if (is_anchor_forward) {
    const int id = current_forward_anchor_type_preference_ - 1;
    ++forward_anchor_type_count_[id];
    bool update_preference = (id != 0 ) &&
      forward_anchor_type_count_[id] > forward_anchor_type_count_[id - 1];
    
    if (update_preference) {
      int temp = forward_anchor_type_count_[id];
      forward_anchor_type_count_[id] = forward_anchor_type_count_[id - 1];
      forward_anchor_type_count_[id - 1] = temp;
      temp = forward_anchor_type_preference_[id];
      forward_anchor_type_preference_[id] = forward_anchor_type_preference_[id - 1];
      forward_anchor_type_preference_[id - 1] = temp;
    }

  } else {
    const int id = current_reverse_anchor_type_preference_ - 1;
    ++reverse_anchor_type_count_[id];
    bool update_preference = (id != 0 ) &&
      reverse_anchor_type_count_[id] > reverse_anchor_type_count_[id - 1];
    
    if (update_preference) {
      int temp = reverse_anchor_type_count_[id];
      reverse_anchor_type_count_[id] = reverse_anchor_type_count_[id - 1];
      reverse_anchor_type_count_[id - 1] = temp;
      temp = reverse_anchor_type_preference_[id];
      reverse_anchor_type_preference_[id] = reverse_anchor_type_preference_[id - 1];
      reverse_anchor_type_preference_[id - 1] = temp;
    } // if
  } // if-else

  RewindRegionTypeList();

  return true;

}

bool SearchRegionType::GetNextRegionType(const bool is_anchor_forward, 
    RegionType* region_type) {

  if (is_anchor_forward) {
    int preference_id = current_forward_anchor_type_preference_;
    if (preference_id == types_) {
      has_gotten_forward_type_ = false; // make SetCurrentTypeSuccess fail
      return false;
    } else {
      // output current type and increase preference_id
      int type_id = forward_anchor_type_preference_[preference_id];
      *region_type = forward_anchor_type_vector_[type_id];
      ++current_forward_anchor_type_preference_;
      has_gotten_forward_type_ = true;
    }
  
  } else {
    int preference_id = current_reverse_anchor_type_preference_;
    if (preference_id == types_) {
      has_gotten_reverse_type_ = false; // make SetCurrentTypeSuccess fail
      return false;
    } else {
      // output current type and increase preference_id
      int type_id = reverse_anchor_type_preference_[preference_id];
      *region_type = reverse_anchor_type_vector_[type_id];
      ++current_reverse_anchor_type_preference_;
      has_gotten_reverse_type_ = true;
    }
  }

  return true;
}

void SearchRegionType::ResetRegionTypeList(void) {
  // reset the preference as defaults
  for (int i = 0; i < types_; ++i) {
    forward_anchor_type_preference_[i] = i;
    reverse_anchor_type_preference_[i] = i;
    forward_anchor_type_count_[i] = 0;
    reverse_anchor_type_count_[i] = 0;
  }
  RewindRegionTypeList();
}
} // namespace
