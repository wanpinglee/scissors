#include "aligner.h"

#include <list>

#include "dataStructures/target_region.h"
#include "dataStructures/optional_tag.h"
#include "utilities/bam/bam_constant.h"
#include "utilities/bam/bam_utilities.h"
#include "utilities/miscellaneous/alignment_filter.h"
#include "utilities/miscellaneous/hashes_collection.h"

namespace Scissors {
namespace {
const uint32_t kMediumSizedIndelMax = 200;
const uint32_t kMediumSizedIndelMin = 20;

void SetTargetSequence(const SearchRegionType::RegionType& region_type, 
                       SR_QueryRegion* query_region) {
  if (region_type.sequence_inverse && region_type.sequence_complement) {
      SR_QueryRegionChangeSeq(query_region, SR_REVERSE_COMP);
    SR_SetStrand(query_region->pOrphan, SR_REVERSE_COMP);
    query_region->isOrphanInversed = FALSE;

  } else if (region_type.sequence_inverse && !region_type.sequence_complement) {
    SR_QueryRegionChangeSeq(query_region, SR_INVERSE);
      
    SR_SetStrand(query_region->pOrphan, SR_INVERSE);
    query_region->isOrphanInversed = TRUE;
 
  } else if (!region_type.sequence_inverse && region_type.sequence_complement) {
    SR_QueryRegionChangeSeq(query_region, SR_COMP);

    SR_SetStrand(query_region->pOrphan, SR_REVERSE_COMP);
    query_region->isOrphanInversed = TRUE;
  
  } else {
    SR_QueryRegionChangeSeq(query_region, SR_FORWARD);
    SR_SetStrand(query_region->pOrphan, SR_FORWARD);
    query_region->isOrphanInversed = FALSE;
  } // end if-else
  
}

void AdjustBamFlag(bam1_t* al_bam_anchor, bam1_t* partial_al1, bam1_t* partial_al2) {
  namespace Constant = BamFlagConstant;
  bool anchor_mate1 = al_bam_anchor->core.flag & Constant::kBamFMate1;
  bool anchor_reverse = al_bam_anchor->core.flag & Constant::kBamFReverse;

  // set anchor's flag
  al_bam_anchor->core.flag &= 0xfff7; // mate is not unmapped
  if (anchor_reverse) al_bam_anchor->core.flag &= 0xffdf; // mate is not reverse
  else al_bam_anchor->core.flag |= Constant::kBamFReverseMate;

  //set partial alignments flags
  int16_t partial_flag = 0;
  partial_flag |= Constant::kBamFPaired;
  if (anchor_reverse) partial_flag |= Constant::kBamFReverseMate;
  else partial_flag |= Constant::kBamFReverse;
  if (anchor_mate1) partial_flag |= Constant::kBamFMate2;
  else partial_flag |= Constant::kBamFMate1;
  partial_al1->core.flag = partial_flag;
  partial_al2->core.flag = partial_flag;
}

bool ExistMediumIndel(const StripedSmithWaterman::Alignment& ssw_al) {
  for (unsigned int i = 0; i < ssw_al.cigar.size(); ++i) {
    uint8_t  op = ssw_al.cigar[i] & 0x0f;
    uint32_t length = 0;
    switch(op) {
      case 1: //I
        length = ssw_al.cigar[i] >> 4;
	if ((kMediumSizedIndelMin <= length) && (length <= kMediumSizedIndelMax))
	  return true;
	break;
      case 2: //D
        length = ssw_al.cigar[i] >> 4;
	if ((kMediumSizedIndelMin <= length) && (length <= kMediumSizedIndelMax))
	  return true;
        break;
      default:
        break;
    } // end switch
  } // end for

  return false;
}

} // unnamed namespace

Aligner::Aligner(const SR_Reference* reference, 
                 const SR_InHashTable* hash_table,
		 const SR_Reference* reference_special,
		 const SR_InHashTable* hash_table_special,
		 const SR_RefHeader* reference_header,
		 const int& fragment_length)
    : search_region_type_()
    , anchor_region_()
    , reference_(reference)
    , hash_table_(hash_table)
    , reference_special_(reference_special)
    , hash_table_special_(hash_table_special)
    , reference_header_(reference_header)
    , banded_sw_aligner_()
    , stripe_sw_aligner_(){
  
  query_region_     = SR_QueryRegionAlloc();
  hashes_           = HashRegionTableAlloc();
  hashes_special_   = HashRegionTableAlloc();
  special_ref_view_ = SR_RefViewAlloc();

  stripe_sw_aligner_.SetGapPenalty(3, 0);

  //hash_length_.fragLen = fragment_length;
  //hash_length_.closeRange = 2000;
  //hash_length_.farRange   = 10000;
}

Aligner::~Aligner() {
  SR_QueryRegionFree(query_region_);
  HashRegionTableFree(hashes_);
  HashRegionTableFree(hashes_special_);
  SR_RefViewFree(special_ref_view_);
}

void Aligner::LoadRegionType(const bam1_t& anchor) {
  uint32_t* cigar = bam1_cigar(&anchor);
  bool is_new_region = anchor_region_.IsNewRegion(cigar,
                                                  anchor.core.n_cigar, 
                                                  anchor.core.pos);
  
  if (is_new_region)
    search_region_type_.ResetRegionTypeList();
  else
    search_region_type_.RewindRegionTypeList();
}

void Aligner::AlignCandidate(const TargetEvent& target_event,
                             const TargetRegion& target_region,
			     const AlignmentFilter& alignment_filter,
			     SR_BamInStreamIter* al_ite,
                             vector<bam1_t*>* alignments) {
  // Since stripe-smith-waterman is used for local search, 
  // closeRange may not be necessary.
  hash_length_.fragLen    = target_region.fragment_length;
  hash_length_.closeRange = target_region.local_window_size;
  hash_length_.farRange   = target_region.discovery_window_size;
  while (SR_QueryRegionLoadPair(query_region_, al_ite) == SR_OK) {
    // TODO@WP: it may be removed later
    if (query_region_->algnType != SR_UNIQUE_ORPHAN) continue;

    Align(target_event, target_region, alignment_filter, query_region_, alignments);
  } // end while

  al_ite = NULL;
}

void Aligner::AlignCandidate(const TargetEvent& target_event,
                             const TargetRegion& target_region,
		             const AlignmentFilter& alignment_filter,
		             const bam1_t& anchor,
		             const bam1_t& target,
		             vector<bam1_t*>* alignments) {
  // Since stripe-smith-waterman is used for local search, 
  // closeRange may not be necessary.
  hash_length_.fragLen    = target_region.fragment_length;
  hash_length_.closeRange = target_region.local_window_size;
  hash_length_.farRange   = target_region.discovery_window_size;
  
  query_region_->pAnchor = (bam1_t*) &anchor;
  query_region_->pOrphan = (bam1_t*) &target;
  Align(target_event, target_region, alignment_filter, query_region_, alignments);
}

void Aligner::Align(const TargetEvent& target_event,
                    const TargetRegion& target_region,
                    const AlignmentFilter& alignment_filter,
		    const SR_QueryRegion* query_region,
		    vector<bam1_t*>* alignments) {
      
  query_region_ = (SR_QueryRegion*) query_region;
  // Convert 4-bit representive sequence into chars
  SR_QueryRegionLoadSeq(query_region_);
  
  // Search the local region
  StripedSmithWaterman::Alignment ssw_al;
  SearchLocalRegion(target_region, &ssw_al);

  //if (target_event.medium_sized_indel) {
    
   // bool medium_indel_found = SearchMediumIndel(alignment_filter, );
  //} else {
    // nothing
  //}

  // store the anchor in output bam
  /*
  bam1_t* al_bam_anchor;
  al_bam_anchor = bam_init1(); // Thread.cpp will free it
  al_bam_anchor = bam_copy1(al_bam_anchor, query_region_->pAnchor);
  alignments->push_back(al_bam_anchor);
  */
 
  if (target_event.special_insertion) {
    Alignment local_al, special_al;
    bool special_found = SearchSpecialReference(alignment_filter, &local_al, &special_al);
    if (special_found) {
	  // store alignments
	  bam1_t *al1_bam, *al2_bam;
	  al1_bam = bam_init1(); // Thread.cpp will free it
	  al2_bam = bam_init1(); // Thread.cpp will free it
	  BamUtilities::ConvertAlignmentToBam1(local_al, *query_region_->pOrphan, al1_bam);
	  BamUtilities::ConvertAlignmentToBam1(special_al, *query_region_->pOrphan, al2_bam);
	  // add optional tags
	  std::list<bam1_t*> al_list;
	  al_list.push_back(al1_bam);
	  al_list.push_back(al2_bam);
	  OptionalTag::AddOptionalTags(*(query_region_->pAnchor), al_list);

	  uint32_t s_pos;
	  int32_t s_ref_id;
	  SR_GetRefFromSpecialPos(special_ref_view_, &s_ref_id, &s_pos, reference_header_, reference_special_, al2_bam->core.pos);
	  al2_bam->core.pos = s_pos;
	  al2_bam->core.tid = s_ref_id;

	  alignments->push_back(al1_bam);
	  alignments->push_back(al2_bam);
	  //AdjustBamFlag(al_bam_anchor, al1_bam, al2_bam);
    } else { // !special_found
      // nothing
    }
  }

      // For normal detection
      /*
      const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
      SearchRegionType::RegionType region_type;
      LoadRegionType(*(query_region_->pAnchor));
      while (search_region_type_.GetNextRegionType(is_anchor_forward, 
                                                   &region_type)) {
        
	// Reverse or complement the sequence if necesary
        SetTargetSequence(region_type, query_region_);
        HashRegionTableInit(hashes_, read_length);
        SR_QueryRegionSetRange(query_region_, &hash_length_, reference_->seqLen,
                               region_type.upstream ? SR_DOWNSTREAM : SR_UPSTREAM);
	HashRegionTableLoad(hashes_, hash_table_, query_region_);

        HashesCollection hashes_collection;
        HashesCollection hashes_collection_special;
	
	hashes_collection.Init(*(hashes_->pBestCloseRegions));
	//printf("\nBefore sorting\n");
	//hashes_collection.Print();
	hashes_collection.SortByLength();
	//printf("\nAfter sorting\n");
	//hashes_collection.Print();
	unsigned int best1, best2;
	bool best_pair_found = hashes_collection.GetBestCoverPair(&best1, &best2);
	Alignment al1, al2;
	if (best_pair_found) {
	  GetAlignment(hashes_collection, best1, &al1);
	  GetAlignment(hashes_collection, best2, &al2);
	  al1.is_seq_inverse    = region_type.sequence_inverse;
	  al2.is_seq_inverse    = region_type.sequence_inverse;
	  al1.is_seq_complement = region_type.sequence_complement;
	  al2.is_seq_complement = region_type.sequence_complement;

	  // store alignments
	  bam1_t *al1_bam, *al2_bam;
	  al1_bam = bam_init1(); // Thread.cpp will free it
	  al2_bam = bam_init1(); // Thread.cpp will free it
	  BamUtilities::ConvertAlignmentToBam1(al1, *query_region_->pOrphan, al1_bam);
	  BamUtilities::ConvertAlignmentToBam1(al2, *query_region_->pOrphan, al2_bam);
	  alignments->push_back(al1_bam);
	  alignments->push_back(al2_bam);
	} else {
	  // nothing
	}
      } // end while
      */
    //} // end while
}

// @function: Gets alignment from the hashes_collection by given id 
//            in hashes_collection
bool Aligner::GetAlignment(
    const HashesCollection& hashes_collection, 
    const unsigned int& id,
    const bool& special,
    const int& read_length,
    const char* read_seq,
    Alignment* al){
  
  if (static_cast<int>(id) >= hashes_collection.GetSize()) {
    return false;
  } else {
    int hash_begin = (hashes_collection.Get(id))->refBegins[0];
    int begin, end;
    GetTargetRefRegion(read_length, hash_begin, special, &begin, &end);
    int ref_length = end - begin + 1;
    const char* ref_seq = GetSequence(begin, special);
    BandedSmithWatermanHashRegion hr;
    hr.reference_begin = hash_begin - begin;
    hr.query_begin = (hashes_collection.Get(id))->queryBegin;
    banded_sw_aligner_.Align(*al, ref_seq, ref_length, read_seq, read_length, hr);
    al->reference_begin += begin;

    return true;
  }
}

// @function Given a pivot (hash_begin) and the length that you want to extend,
//           the function will set the valid begin and end after extending.
// @param  extend_length The length that you want to extend.
// @param  hash_begin    A pivot that you want to extend the region according to.
// @param  special       In the special reference?
// @param  begin         The resultant begin
// @param  end           The resultant end
void Aligner::GetTargetRefRegion(const int& extend_length, const int& hash_begin, 
    const bool& special, int* begin, int* end) {
  uint32_t pos = hash_begin;
  int32_t ref_id = 0;
  if (special)
    SR_GetRefFromSpecialPos(special_ref_view_, &ref_id, &pos, reference_header_, reference_special_, hash_begin);

  int seq_length = special ? special_ref_view_->seqLen : reference_->seqLen;

  int forward_shift;
  if (pos < static_cast<uint32_t>(extend_length)) forward_shift = pos;
  else forward_shift = extend_length;

  int backward_shift;
  if ((pos + extend_length + 1) > seq_length) backward_shift = seq_length - pos - 1;
  else backward_shift = extend_length;
  
  *begin = hash_begin - forward_shift;
  *end   = hash_begin + backward_shift;
}

bool Aligner::SearchSpecialReference(const AlignmentFilter& alignment_filter,
                                     Alignment* local_al, Alignment* special_al)
{

  SearchRegionType::RegionType region_type;
  const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
  // Loads hashes for the first partial alignment
  search_region_type_.GetStandardType(is_anchor_forward, &region_type);
  SetTargetSequence(region_type, query_region_);
  int read_length = query_region_->pOrphan->core.l_qseq;
  HashRegionTableInit(hashes_, read_length);
  SR_QueryRegionSetRange(query_region_, &hash_length_, reference_->seqLen,
                         region_type.upstream ? SR_DOWNSTREAM : SR_UPSTREAM);
  HashRegionTableLoad(hashes_, hash_table_, query_region_);
  HashesCollection hashes_collection;
  hashes_collection.Init(*(hashes_->pBestCloseRegions));
     
  // Loads special hashes
  HashRegionTableInit(hashes_special_, read_length);
  SR_QueryRegionSetRangeSpecial(query_region_, reference_special_->seqLen);
  HashRegionTableLoad(hashes_special_, hash_table_special_, query_region_);
  HashesCollection hashes_collection_special;
  hashes_collection_special.Init(*(hashes_special_->pBestCloseRegions));
  //hashes_collection_special.Print();

  // Gets the best pair of hashes
  unsigned int best1, best2;
  bool best_pair_found = false;
  if (is_anchor_forward)
    best_pair_found = hashes_collection.GetBestCoverPair(&hashes_collection_special, &best1, &best2);
  else 
    best_pair_found = hashes_collection_special.GetBestCoverPair(&hashes_collection, &best2, &best1);
      
  // If a pair of hashes is found
  if (best_pair_found) {
    const char* read_seq = query_region_->orphanSeq;
    GetAlignment(hashes_collection, best1, false, read_length, read_seq, local_al); // non-special
    GetAlignment(hashes_collection_special, best2, true, read_length, read_seq, special_al); // special
	
    namespace filter_app = AlignmentFilterApplication;
    bool trimming_al_okay = true;
    //trimming_al_okay &= filter_app::TrimAlignment(alignment_filter, &al1);
    //trimming_al_okay &= filter_app::TrimAlignment(alignment_filter, &al2);
    filter_app::TrimAlignment(alignment_filter, local_al);
    filter_app::TrimAlignment(alignment_filter, special_al);
    trimming_al_okay &= ((local_al->reference.size() > 0) && (special_al->reference.size() > 0));

    bool passing_filter = true;
    passing_filter &= filter_app::FilterByMismatch(alignment_filter, *local_al);
    passing_filter &= filter_app::FilterByMismatch(alignment_filter, *special_al);
    passing_filter &= filter_app::FilterByAlignedBaseThreshold(alignment_filter, *local_al, read_length);
    passing_filter &= filter_app::FilterByAlignedBaseThreshold(alignment_filter, *special_al, read_length);

    if (trimming_al_okay && passing_filter) {
      local_al->  is_seq_inverse    = region_type.sequence_inverse;
      special_al->is_seq_inverse    = region_type.sequence_inverse;
      local_al->  is_seq_complement = region_type.sequence_complement;
      special_al->is_seq_complement = region_type.sequence_complement;
      return true;
    }
    else {
      return false;
    }
  } else {
    // !best_pair_found
    return false;
  }
}

void Aligner::SearchLocalRegion(const TargetRegion& target_region,
                                StripedSmithWaterman::Alignment* ssw_al) {
  // Get the standard region type
  SearchRegionType::RegionType region_type;
  const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
  search_region_type_.GetStandardType(is_anchor_forward, &region_type);

  // Calculate the pivot position
  int anchor_pos = query_region_->pAnchor->core.pos;
  int pivot = region_type.upstream ? anchor_pos + target_region.fragment_length
                                   : anchor_pos - target_region.fragment_length;
  if (pivot < 0) pivot = 0;

  // Get the region according to the pivot and local_window_size
  int begin, end;
  bool special = false;
  GetTargetRefRegion(target_region.local_window_size, pivot, special, &begin, &end);

  // Get the read sequence
  SetTargetSequence(region_type, query_region_);
  string read_seq;
  read_seq.assign(query_region_->orphanSeq, query_region_->pOrphan->core.l_qseq);

  // Apply SSW to the region
  int ref_length = end - begin + 1;
  const char* ref_seq = GetSequence(begin, special);
  StripedSmithWaterman::Filter filter;
  filter.distance_filter = kMediumSizedIndelMax;
  stripe_sw_aligner_.Align(read_seq.c_str(), ref_seq, ref_length, filter, ssw_al);

  /*
  if (!ssw_al.cigar_string.empty()) {
    fprintf(stderr, "%u\t%u\t%s\t%s\t%u\t%u\t%u\t%u\n", query_region_->pAnchor->core.tid, anchor_pos, read_seq.c_str(), ssw_al.cigar_string.c_str(), ssw_al.query_begin, ssw_al.query_end, ssw_al.ref_begin, ssw_al.ref_end);
  }
  */
}

inline const char* Aligner::GetSequence(const size_t& start, const bool& special) const {
  if (special)
    return (reference_special_->sequence + start);
  else
    return (reference_->sequence + start);
}
} //namespace
