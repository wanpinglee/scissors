#include "aligner.h"

#include <list>

#include "dataStructures/target_event.h"
#include "dataStructures/target_region.h"
#include "dataStructures/optional_tag.h"
#include "utilities/bam/bam_constant.h"
#include "utilities/bam/bam_utilities.h"
#include "utilities/miscellaneous/alignment_filter.h"
#include "utilities/miscellaneous/alignment_collection.h"
#include "utilities/miscellaneous/hashes_collection.h"

namespace Scissors {
namespace {

const Scissors::TargetEvent kMediumIndel(false, true, false);
const Scissors::TargetEvent kSpecialInsertion(true, false, false);
const Scissors::TargetEvent kInsertion(false, false, true);


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

void AdjustBamFlag(const bam1_t& al_bam_anchor, bam1_t* partial_al1) {
  namespace Constant = BamFlagConstant;
  const bool anchor_mate1 = al_bam_anchor.core.flag & Constant::kBamFMate1;
  bool anchor_reverse = al_bam_anchor.core.flag & Constant::kBamFReverse;

  // set anchor's flag
  //al_bam_anchor->core.flag &= 0xfff7; // mate is not unmapped
  //if (anchor_reverse) al_bam_anchor->core.flag &= 0xffdf; // mate is not reverse
  //else al_bam_anchor->core.flag |= Constant::kBamFReverseMate;

  //set partial alignments flags
  int16_t partial_flag = 0;
  partial_flag |= Constant::kBamFPaired;
  if (anchor_reverse) partial_flag |= Constant::kBamFReverseMate;
  else partial_flag |= Constant::kBamFReverse;
  if (anchor_mate1) partial_flag |= Constant::kBamFMate2;
  else partial_flag |= Constant::kBamFMate1;
  partial_al1->core.flag = partial_flag;
}

bool FindMediumIndel(const StripedSmithWaterman::Alignment& ssw_al,
                     unsigned int* indel_cigar_id) {
  uint32_t longest = 0;
  bool found = false;
  for (unsigned int i = 0; i < ssw_al.cigar.size(); ++i) {
    uint8_t  op = ssw_al.cigar[i] & 0x0f;
    uint32_t length = 0;
    switch(op) {
      case 1: // I
      case 2: // D
        length = ssw_al.cigar[i] >> 4;
	if (length > longest) {
	  longest = length;
	  *indel_cigar_id = i;
	}
	found = true;
	break;
      default:
        break;
    } // end switch
  } // end for

  return found;
}

bool DetermineMediumIndelBreakpoint(
    const StripedSmithWaterman::Alignment& ssw_al,
    const unsigned int& indel_cigar_id,
    const AlignmentFilter& alignment_filter,
    const int& read_length) {
  // The statement should not be true.
  if (ssw_al.cigar.empty()) return false;

  bool insertion = ((ssw_al.cigar[indel_cigar_id] & 0x0f) == 1) ? true : false;

  unsigned int aligned_base_threshold = floor(read_length * alignment_filter.aligned_base_rate);
  bool pass_filter = false;
  for(unsigned int i = 0; i < indel_cigar_id; ++i) {
    uint8_t  op = ssw_al.cigar[i] & 0x0f;
    uint32_t length = ssw_al.cigar[i] >> 4;
    // Cigar operation: M
    if ((op == 0) && (length > aligned_base_threshold)) {
      pass_filter = true;
      break;
    }
  }

  // The first partial doesn't exist to support the detected deletion.
  // So, return false. For deletions, both side partital alignments should flank.
  if (!pass_filter && !insertion) {
    return false;
  }

  // The first partial is found to support the detected insertion.
  // So, return true;
  if (pass_filter && insertion) return true;

  pass_filter = false;
  for(unsigned int i = indel_cigar_id + 1; i < ssw_al.cigar.size(); ++i) {
    uint8_t  op = ssw_al.cigar[i] & 0x0f;
    uint32_t length = ssw_al.cigar[i] >> 4;
    //fprintf(stderr, "%u\t%u\t%u\t", i, op, length);
    if ((op == 0) && (length > aligned_base_threshold)) {
      pass_filter = true;
      break;
    }
  }
  //fprintf(stderr, "%c\n", pass_filter?'T':'F');
  return pass_filter;
}

bool PassMismatchFilter(
    const StripedSmithWaterman::Alignment& ssw_al,
    const AlignmentFilter& alignment_filter,
    const int& event_length) {
  int bases = 0;
  for (unsigned int i = 0; i < ssw_al.cigar.size(); ++i) {
    uint8_t  op = ssw_al.cigar[i] & 0x0f;
    uint32_t length = 0;
    switch(op) {
      case 0: // M
      case 1: // I
        length = ssw_al.cigar[i] >> 4;
	bases += length;
	break;
      default:
        break;
    } // end switch
  } // end for

  float allowed_mismatches = ceil(bases * alignment_filter.allowed_mismatch_rate);
  int mismatches = ssw_al.mismatches - event_length;
  if (mismatches < static_cast<int>(allowed_mismatches)) return true;
  else return false;
}


void StoreAlignment(
    const TargetEvent& event,
    const vector <StripedSmithWaterman::Alignment*>& ssw_als,
    const vector <Alignment*>& common_als,
    const bam1_t& anchor,
    const bam1_t& target,
    vector<bam1_t*>* alignments) {
  // Stores alignments that are obtained by SSW
  for (vector <StripedSmithWaterman::Alignment*>::const_iterator ite = ssw_als.begin();
      ite != ssw_als.end(); ++ite) {
    bam1_t *al_bam;
    al_bam = bam_init1(); // Thread.cpp will free it
    //const bool is_anchor_forward = !bam1_strand(&anchor);
    //fprintf(stderr,"%c%c\n", (*ite)->is_reverse?'T':'F', (*ite)->is_complement?'T':'F');
    BamUtilities::ConvertAlignmentToBam1(**ite, target, (*ite)->is_reverse, (*ite)->is_complement, al_bam);
    AdjustBamFlag(anchor, al_bam);
    OptionalTag::AddOptionalTags(anchor, al_bam);
    alignments->push_back(al_bam);
  }

  for (vector <Alignment*>::const_iterator ite = common_als.begin();
      ite != common_als.end(); ++ite) {
    bam1_t *al_bam;
    al_bam = bam_init1(); // Thread.cpp will free it
    BamUtilities::ConvertAlignmentToBam1(**ite, target, al_bam);
    AdjustBamFlag(anchor, al_bam);
    OptionalTag::AddOptionalTags(anchor, al_bam);
    alignments->push_back(al_bam);
  }
}

bool CheckSetting(const SR_Reference* reference, const Technology technology) {
  if (reference == NULL) return false;
  if (technology == TECH_NONE) return false;
  return true;
}

} // unnamed namespace

Aligner::Aligner() 
    : search_region_type_()
    , anchor_region_()
    , technology_(TECH_NONE)
    , reference_(NULL)
    , hash_table_(NULL)
    , reference_special_(NULL)
    , hash_table_special_(NULL)
    , reference_header_(NULL)
    , query_region_(NULL)
    , hashes_(NULL)
    , hashes_special_(NULL)
    , hash_length_()
    , special_ref_view_()
    , banded_sw_aligner_()
    , stripe_sw_indel_()
    , stripe_sw_normal_() {
  query_region_     = SR_QueryRegionAlloc();
  hashes_           = HashRegionTableAlloc();
  hashes_special_   = HashRegionTableAlloc();
  special_ref_view_ = SR_RefViewAlloc();

  stripe_sw_indel_.Clear();
  stripe_sw_indel_.ReBuild(10,20,20,1);
}

Aligner::Aligner(const SR_Reference* reference, 
                 const SR_InHashTable* hash_table,
		 const SR_Reference* reference_special,
		 const SR_InHashTable* hash_table_special,
		 const SR_RefHeader* reference_header,
		 const Technology& technology)
    : search_region_type_()
    , anchor_region_()
    , technology_(technology)
    , reference_(reference)
    , hash_table_(hash_table)
    , reference_special_(reference_special)
    , hash_table_special_(hash_table_special)
    , reference_header_(reference_header)
    , query_region_(NULL)
    , hashes_(NULL)
    , hashes_special_(NULL)
    , hash_length_()
    , special_ref_view_()
    , banded_sw_aligner_()
    , stripe_sw_indel_()
    , stripe_sw_normal_() {
  
  query_region_     = SR_QueryRegionAlloc();
  hashes_           = HashRegionTableAlloc();
  hashes_special_   = HashRegionTableAlloc();
  special_ref_view_ = SR_RefViewAlloc();

  stripe_sw_indel_.Clear();
  stripe_sw_indel_.ReBuild(10,20,20,1);
  //stripe_sw_aligner_.SetGapPenalty(4, 0);

  //hash_length_.fragLen = fragment_length;
  //hash_length_.closeRange = 2000;
  //hash_length_.farRange   = 10000;
}

bool Aligner::SetReference(
    const SR_Reference*   reference,
    const SR_InHashTable* hash_table,
    const Technology&     technology,
    const SR_Reference*   reference_special,
    const SR_InHashTable* hash_table_special,
    const SR_RefHeader*   reference_header) {
  reference_          = reference;
  hash_table_         = hash_table;
  technology_         = technology;
  reference_special_  = reference_special;
  hash_table_special_ = hash_table_special;
  reference_header_   = reference_header;
  return true;
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

bool Aligner::AlignCandidate(const TargetEvent& target_event,
                             const TargetRegion& target_region,
			     const AlignmentFilter& alignment_filter,
			     SR_BamInStreamIter* al_ite,
                             vector<bam1_t*>* alignments) {
  if (!CheckSetting(reference_, technology_)) {
    while (SR_QueryRegionLoadPair(query_region_, al_ite) == SR_OK); // empty the buffer
    al_ite = NULL;
    return false;
  }
  
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
  return true;
}

bool Aligner::AlignCandidate(const TargetEvent& target_event,
                             const TargetRegion& target_region,
		             const AlignmentFilter& alignment_filter,
		             const bam1_t& anchor,
		             const bam1_t& target,
		             vector<bam1_t*>* alignments) {
  if (!CheckSetting(reference_, technology_)) {
    return false;
  }

  // Since stripe-smith-waterman is used for local search, 
  // closeRange may not be necessary.
  hash_length_.fragLen    = target_region.fragment_length;
  hash_length_.closeRange = target_region.local_window_size;
  hash_length_.farRange   = target_region.discovery_window_size;
  
  query_region_->pAnchor = (bam1_t*) &anchor;
  query_region_->pOrphan = (bam1_t*) &target;
  Align(target_event, target_region, alignment_filter, query_region_, alignments);

  return true;
}

void Aligner::Align(const TargetEvent& target_event,
                    const TargetRegion& target_region,
                    const AlignmentFilter& alignment_filter,
		    const SR_QueryRegion* query_region,
		    vector<bam1_t*>* alignments) {
      
  query_region_ = (SR_QueryRegion*) query_region;
  // Convert 4-bit representive sequence into chars
  SR_QueryRegionLoadSeq(query_region_);
  
  AlignmentCollection al_collection;
  StripedSmithWaterman::Alignment indel_al;
  // Try to align for medium-sized INDELs
  if (target_event.medium_sized_indel) {
    bool medium_indel_found = 
        SearchMediumIndel(target_region, alignment_filter, &indel_al);
    if (medium_indel_found) { // push the event and its corresponding alignments in the collection
      al_collection.PushANewEvent(kMediumIndel);
      al_collection.PushAlignment(indel_al);
    } else { // !medium_indel_found
      // nothing
    }
  } else { // !target_event.medium_sized_indel
    // nothing
  }

  // Try to align for special insertions
  StripedSmithWaterman::Alignment local_al;
  Alignment special_al;
  if (target_event.special_insertion) {
    bool special_found = SearchSpecialReference(target_region, alignment_filter, &special_al);
    if (special_found) {
      //local_al.reference_id = query_region_->pOrphan->core.tid;
      // adjust reference id for alignments in special insertions
      uint32_t s_pos;
      uint32_t s_length = special_al.reference_end - special_al.reference_begin;
      int32_t s_ref_id;
      SR_GetRefFromSpecialPos(special_ref_view_, &s_ref_id, &s_pos, reference_header_,
          reference_special_, special_al.reference_begin);
      special_al.reference_begin = s_pos;
      special_al.reference_end = s_pos + s_length;
      special_al.reference_id = s_ref_id;
	  
      // push the event and its corresponding alignments in the collection
      al_collection.PushANewEvent(kSpecialInsertion);
      al_collection.PushAlignment(local_al);
      al_collection.PushAlignment(special_al);
    } else { // !special_found
      // nothing
    }
  }

  // Pick up the most supported event and its corresponding alignments
  TargetEvent best_event;
  vector <StripedSmithWaterman::Alignment*> ssw_al_for_best_event;
  vector <Alignment*> common_al_for_best_event;
  al_collection.GetMostConfidentEvent(&best_event, &ssw_al_for_best_event, &common_al_for_best_event);
  StoreAlignment(best_event, ssw_al_for_best_event, common_al_for_best_event, 
      *query_region_->pAnchor, *query_region_->pOrphan, alignments);
}

inline bool Aligner::ResetIndelSmithWatermanScore(
    const uint8_t& match_score,
    const uint8_t& mismatch_penalty,
    const uint8_t& gap_opening_penalty,
    const uint8_t& gap_extending_penalty) {
  stripe_sw_indel_.Clear();
  return stripe_sw_indel_.ReBuild(10,20,20,1);
}

inline bool Aligner::ResetLocalSmithWatermanScore(
    const uint8_t& match_score,
    const uint8_t& mismatch_penalty,
    const uint8_t& gap_opening_penalty,
    const uint8_t& gap_extending_penalty) {
  stripe_sw_normal_.Clear();
  return stripe_sw_normal_.ReBuild(10,20,20,1);
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

bool Aligner::SearchLocalPartial(const TargetRegion& target_region,
                                 const AlignmentFilter& alignment_filter,
                                 StripedSmithWaterman::Alignment* local_al) {
  namespace Constant = BamFlagConstant;
  const bool is_anchor_mate1 = query_region_->pAnchor->core.flag & Constant::kBamFMate1;
  const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
  search_region_type_.SetTechnologyAndAnchorMate1(technology_, is_anchor_mate1);
  SearchRegionType::RegionType region_type;

  // ============================================
  // Loads hashes for the first partial alignment
  // ============================================
  search_region_type_.GetStandardType(is_anchor_forward, &region_type);
  SetTargetSequence(region_type, query_region_);

  // Calculate the pivot position and target region
  int anchor_pos = query_region_->pAnchor->core.pos;
  int pivot = region_type.upstream ? anchor_pos + target_region.fragment_length
                                   : anchor_pos - target_region.fragment_length;
  if (pivot < 0) pivot = 0;

  // Get the region according to the pivot and local_window_size
  int begin, end;
  bool special = false;
  GetTargetRefRegion(target_region.local_window_size, pivot, special, &begin, &end);
  if (end < begin) return false;

  // Get the read sequence
  SetTargetSequence(region_type, query_region_);
  string read_seq;
  read_seq.assign(query_region_->orphanSeq, query_region_->pOrphan->core.l_qseq);

  // Apply SSW to the region
  int ref_length = end - begin + 1;
  const char* ref_seq = GetSequence(begin, special);
  StripedSmithWaterman::Filter filter;
  // If the difference between beginnng and ending is larger than distance_filter,
  // then we don't think that it's medium-sized indels and don't need to trace
  // the alignment.
  filter.distance_filter = read_seq.size();
  stripe_sw_normal_.Align(read_seq.c_str(), ref_seq, ref_length, filter, local_al);
  local_al->is_reverse = region_type.sequence_inverse;
  local_al->is_complement = region_type.sequence_complement;
  // Return false since no alignment is found
  //fprintf(stderr, "%u\t%u\t%s\n", indel_al->ref_end, indel_al->ref_begin, indel_al->cigar_string.c_str());
  if (local_al->cigar.empty()) return false;
  namespace filter_app = AlignmentFilterApplication;
  filter_app::TrimAlignment(alignment_filter, local_al);
  if (local_al->cigar.size() == 0) return false;
  else {
    if (!PassMismatchFilter(*local_al, alignment_filter, 0)) return false;
  }

  // Adjust the positions
  local_al->ref_begin += begin;
  local_al->ref_end   += begin;
  local_al->ref_end_next_best += begin;

  return true;
}

bool Aligner::SearchSpecialReference(const TargetRegion& target_region,
                                     const AlignmentFilter& alignment_filter,
				     Alignment* special_al)
{
  
  namespace Constant = BamFlagConstant;
  const bool is_anchor_mate1 = query_region_->pAnchor->core.flag & Constant::kBamFMate1;
  const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
  search_region_type_.SetTechnologyAndAnchorMate1(technology_, is_anchor_mate1);
  SearchRegionType::RegionType region_type;
  search_region_type_.GetStandardType(is_anchor_forward, &region_type);
  SetTargetSequence(region_type, query_region_);
  int read_length = query_region_->pOrphan->core.l_qseq;

  // Get the read sequence
  SetTargetSequence(region_type, query_region_);
  string read_seq;
  read_seq.assign(query_region_->orphanSeq, query_region_->pOrphan->core.l_qseq);

  // ====================
  // Loads special hashes
  // ====================
  HashRegionTableInit(hashes_special_, read_length);
  SR_QueryRegionSetRangeSpecial(query_region_, reference_special_->seqLen);
  HashRegionTableLoad(hashes_special_, hash_table_special_, query_region_);
  HashesCollection hashes_collection_special;
  hashes_collection_special.Init(*(hashes_special_->pBestCloseRegions));
  //hashes_collection_special.Print();
  hashes_collection_special.SortByLength();
  if (hashes_collection_special.Get(hashes_collection_special.GetSize() - 1) == NULL) return false;
  if (hashes_collection_special.Get(hashes_collection_special.GetSize() - 1)->length == 0) return false;

  const int best2 = hashes_collection_special.GetSize() - 1;
  GetAlignment(hashes_collection_special, best2, true, read_length, read_seq.c_str(), special_al);
  bool trimming_al_okay = true;
  namespace filter_app = AlignmentFilterApplication;
  filter_app::TrimAlignment(alignment_filter, special_al);
  trimming_al_okay &= (special_al->reference.size() > 0);

  bool passing_filter = true;
  passing_filter &= filter_app::FilterByMismatch(alignment_filter, *special_al);
  passing_filter &= filter_app::FilterByAlignedBaseThreshold(alignment_filter, *special_al, read_length);
  
  if (trimming_al_okay && passing_filter) {
    special_al->is_seq_inverse    = region_type.sequence_inverse;
    special_al->is_seq_complement = region_type.sequence_complement;
    return true;
  } else {
    return false;
  }
}

bool Aligner::SearchMediumIndel(const TargetRegion& target_region,
                                const AlignmentFilter& alignment_filter,
                                StripedSmithWaterman::Alignment* indel_al) {
  
  namespace Constant = BamFlagConstant;
  const bool is_anchor_mate1 = query_region_->pAnchor->core.flag & Constant::kBamFMate1;
  const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
  search_region_type_.SetTechnologyAndAnchorMate1(technology_, is_anchor_mate1);
  SearchRegionType::RegionType region_type;

  // ============================
  // Get the standard region type
  // ============================
  search_region_type_.GetStandardType(is_anchor_forward, &region_type);

  //fprintf(stderr,"%s\n", bam1_qname(query_region_->pAnchor));
  //fprintf(stderr, "%c%c%c\n", region_type.upstream?'T':'F', region_type.sequence_inverse?'T':'F', region_type.sequence_complement?'T':'F');

  // ==============================================
  // Calculate the pivot position and target region
  // ==============================================
  int anchor_pos = query_region_->pAnchor->core.pos;
  int pivot = region_type.upstream ? anchor_pos + target_region.fragment_length
                                   : anchor_pos - target_region.fragment_length;
  if (pivot < 0) pivot = 0;

  // Get the region according to the pivot and local_window_size
  int begin, end;
  bool special = false;
  GetTargetRefRegion(target_region.local_window_size, pivot, special, &begin, &end);
  //fprintf(stderr, "%u\t%u\n", begin, end);
  if (end < begin) return false;

  // =====================
  // Get the read sequence
  // =====================
  SetTargetSequence(region_type, query_region_);
  string read_seq;
  read_seq.assign(query_region_->orphanSeq, query_region_->pOrphan->core.l_qseq);
  //fprintf(stderr, "%s\n", read_seq.c_str());

  // =======================
  // Apply SSW to the region
  // =======================
  int ref_length = end - begin + 1;
  const char* ref_seq = GetSequence(begin, special);
  //fprintf(stderr, "%d\t%u\n", reference_->id, reference_->seqLen);
  //for (unsigned int i = 0; i < 10; ++i)
  //  fprintf(stderr, "%c", *(ref_seq+i));
  //fprintf(stderr, "\n");
  StripedSmithWaterman::Filter filter;
  // If the difference between beginnng and ending is larger than distance_filter,
  // then we don't think that it's medium-sized indels and don't need to trace
  // the alignment.
  filter.distance_filter = read_seq.size() * 5;
  stripe_sw_indel_.Align(read_seq.c_str(), ref_seq, ref_length, filter, indel_al);
  indel_al->is_reverse = region_type.sequence_inverse;
  indel_al->is_complement = region_type.sequence_complement;
  // Return false since no alignment is found
  //fprintf(stderr, "%u\t%u\t%s\n", indel_al->ref_end, indel_al->ref_begin, indel_al->cigar_string.c_str());
  if (indel_al->cigar.empty()) return false;

  // Adjust the positions
  indel_al->ref_begin += begin;
  indel_al->ref_end   += begin;
  indel_al->ref_end_next_best += begin;

  // ==================================
  // Check if medium-sized INDELs exist
  // ==================================
  unsigned int cigar_id = 0;
  bool indel_found = FindMediumIndel(*indel_al, &cigar_id);
  indel_found &= DetermineMediumIndelBreakpoint(*indel_al, cigar_id, alignment_filter, read_seq.size());
  int event_length = indel_al->cigar[cigar_id] >> 4;

  //return indel_found;

  if (indel_found) {
    if (indel_al->cigar.size() > 7) {
      AlignmentFilterApplication::TrimAlignment(alignment_filter, indel_al);
      if (indel_al->cigar.size() == 0) return false;
      else return PassMismatchFilter(*indel_al, alignment_filter, event_length);
    } else { // deletion
      return true;
    }
    //fprintf(stderr, "%s\t%u\n", indel_al->cigar_string.c_str(), cigar_id);
    //return DetermineMediumIndelBreakpoint(*indel_al, cigar_id, alignment_filter, read_seq.size());
  } else { // !indel_found
    return false;
  }


}

inline const char* Aligner::GetSequence(const size_t& start, const bool& special) const {
  if (special)
    return (reference_special_->sequence + start);
  else
    return (reference_->sequence + start);
}
} //namespace
