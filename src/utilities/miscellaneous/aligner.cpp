#include "aligner.h"
#include "utilities/bam/bam_constant.h"
#include "utilities/bam/bam_utilities.h"
#include "utilities/miscellaneous/alignment_filter.h"
#include "utilities/miscellaneous/hashes_collection.h"

namespace {
void SetTargetSequence(const SearchRegionType::RegionType& region_type, 
                       SR_QueryRegion* query_region) {
  //const bool forward    = !bam1_strand(query_region->pOrphan);
  //const bool inverse    = (forward && (query_region->isOrphanInversed == TRUE))
  //                      ||(!forward && (query_region->isOrphanInversed == FALSE));
  //const bool complement = !forward;
  
  if (region_type.sequence_inverse && region_type.sequence_complement) {
   // if (!inverse && !complement)
      SR_QueryRegionChangeSeq(query_region, SR_REVERSE_COMP);
   // else if (!inverse)
   //   SR_QueryRegionChangeSeq(query_region, SR_INVERSE);
   // else if (!complement)
   //   SR_QueryRegionChangeSeq(query_region, SR_COMP);

    SR_SetStrand(query_region->pOrphan, SR_REVERSE_COMP);
    query_region->isOrphanInversed = FALSE;

  } else if (region_type.sequence_inverse && !region_type.sequence_complement) {
   // if (!inverse && complement)
   //   SR_QueryRegionChangeSeq(query_region, SR_REVERSE_COMP);
   // else if (!inverse)
      SR_QueryRegionChangeSeq(query_region, SR_INVERSE);
   // else if (complement)
   //   SR_QueryRegionChangeSeq(query_region, SR_COMP);
      
    SR_SetStrand(query_region->pOrphan, SR_INVERSE);
    query_region->isOrphanInversed = TRUE;
 
  } else if (!region_type.sequence_inverse && region_type.sequence_complement) {
   // if (inverse && !complement)
   //   SR_QueryRegionChangeSeq(query_region, SR_REVERSE_COMP);
   // else if (inverse)
   //   SR_QueryRegionChangeSeq(query_region, SR_INVERSE);
   // else if (!complement)
      SR_QueryRegionChangeSeq(query_region, SR_COMP);

    SR_SetStrand(query_region->pOrphan, SR_REVERSE_COMP);
    query_region->isOrphanInversed = TRUE;
  
  } else {
   // if (inverse && complement)
   //   SR_QueryRegionChangeSeq(query_region, SR_REVERSE_COMP);
   // else if (inverse)
   //   SR_QueryRegionChangeSeq(query_region, SR_INVERSE);
   // else if (complement)
   //   SR_QueryRegionChangeSeq(query_region, SR_COMP);

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
}

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
    , sw_aligner_(){
  
  query_region_     = SR_QueryRegionAlloc();
  hashes_           = HashRegionTableAlloc();
  hashes_special_   = HashRegionTableAlloc();
  special_ref_view_ = SR_RefViewAlloc();

  hash_length_.fragLen = fragment_length;
  hash_length_.closeRange = 2000;
  hash_length_.farRange   = 10000;
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


void Aligner::AlignCandidate(const bool& detect_special,
                             const AlignmentFilter& alignment_filter,
			     SR_BamInStreamIter* al_ite,
                             vector<bam1_t*>* alignments) {
    while (SR_QueryRegionLoadPair(query_region_, al_ite) == SR_OK) {
      
      // TODO@WP: it may be removed later
      if (query_region_->algnType != SR_UNIQUE_ORPHAN) continue;
      
      const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
      // Convert 4-bit representive sequence into chars
      SR_QueryRegionLoadSeq(query_region_);

      LoadRegionType(*(query_region_->pAnchor));

      SearchRegionType::RegionType region_type;
      HashesCollection hashes_collection;
      HashesCollection hashes_collection_special;
      
      // store the anchor in output bam
      /*
      bam1_t* al_bam_anchor;
      al_bam_anchor = bam_init1(); // Thread.cpp will free it
      al_bam_anchor = bam_copy1(al_bam_anchor, query_region_->pAnchor);
      alignments->push_back(al_bam_anchor);
      */

      // For MEI first
      search_region_type_.GetStandardType(is_anchor_forward, &region_type);
      SetTargetSequence(region_type, query_region_);
      int read_length = query_region_->pOrphan->core.l_qseq;
      HashRegionTableInit(hashes_, read_length);
      SR_QueryRegionSetRange(query_region_, &hash_length_, reference_->seqLen,
                             region_type.upstream ? SR_DOWNSTREAM : SR_UPSTREAM);
      HashRegionTableLoad(hashes_, hash_table_, query_region_);
      hashes_collection.Init(*(hashes_->pBestCloseRegions));
     
      // get special hashes
      HashRegionTableInit(hashes_special_, read_length);
      SR_QueryRegionSetRangeSpecial(query_region_, reference_special_->seqLen);
      HashRegionTableLoad(hashes_special_, hash_table_special_, query_region_);
      hashes_collection_special.Init(*(hashes_special_->pBestCloseRegions));
      //hashes_collection_special.Print();

      // get best cover
      unsigned int best1, best2;
      bool best_pair_found = false;
      if (is_anchor_forward)
        best_pair_found = hashes_collection.GetBestCoverPair(&hashes_collection_special, &best1, &best2);
      else 
        best_pair_found = hashes_collection_special.GetBestCoverPair(&hashes_collection, &best2, &best1);
      
      Alignment al1, al2;
      if (best_pair_found) {
	const char* read_seq = query_region_->orphanSeq;
	GetAlignment(hashes_collection, best1, false, read_length, read_seq, &al1); // non-special
	GetAlignment(hashes_collection_special, best2, true, read_length, read_seq, &al2); // special
	
	namespace filter_app = AlignmentFilterApplication;
	bool trimming_al_okay = true;
	trimming_al_okay &= filter_app::TrimAlignment(alignment_filter, &al1);
	trimming_al_okay &= filter_app::TrimAlignment(alignment_filter, &al2);
	trimming_al_okay &= ((al1.reference.size() > 0) && (al2.reference.size() > 0));

	bool passing_filter = true;
	passing_filter &= filter_app::FilterByMismatch(alignment_filter, al1);
	passing_filter &= filter_app::FilterByMismatch(alignment_filter, al2);
	passing_filter &= filter_app::FilterByAlignedBaseThreshold(alignment_filter, al1, read_length);
	passing_filter &= filter_app::FilterByAlignedBaseThreshold(alignment_filter, al2, read_length);

	if (trimming_al_okay && passing_filter) {
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
	  uint32_t s_pos;
	  int32_t s_ref_id;
	  SR_GetRefFromSpecialPos(special_ref_view_, &s_ref_id, &s_pos, reference_header_, reference_special_, al2_bam->core.pos);
	  al2_bam->core.pos = s_pos;
	  al2_bam->core.tid = s_ref_id;

	  alignments->push_back(al1_bam);
	  alignments->push_back(al2_bam);
	  //AdjustBamFlag(al_bam_anchor, al1_bam, al2_bam);
	} else { // !(trimming_al_okay && passing_filter)
	  // nothing
	}
      } else { // !(best_pair_found)
        // nothing
      }
      
      

      // For normal detection
      /*
      while (search_region_type_.GetNextRegionType(is_anchor_forward, 
                                                   &region_type)) {
        
	// Reverse or complement the sequence if necesary
        SetTargetSequence(region_type, query_region_);
        HashRegionTableInit(hashes_, read_length);
        SR_QueryRegionSetRange(query_region_, &hash_length_, reference_->seqLen,
                               region_type.upstream ? SR_DOWNSTREAM : SR_UPSTREAM);
	HashRegionTableLoad(hashes_, hash_table_, query_region_);

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
    } // end while

    al_ite = NULL;
}

//@description:
//  Gets alignment from the hashes_collection by given id in hashes_collection
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
    sw_aligner_.Align(*al, ref_seq, ref_length, read_seq, read_length, hr);
    al->reference_begin += begin;
    //al->TrimAlignment();
    
    //al->reference_begin = (hashes_collection.Get(id))->refBegins[0];
    //al->reference_end   = (hashes_collection.Get(id))->refBegins[0] + (hashes_collection.Get(id))->length - 1;
    //al->query_begin     = (hashes_collection.Get(id))->queryBegin;
    //al->query_end       = (hashes_collection.Get(id))->queryBegin + (hashes_collection.Get(id))->length - 1;
    //const char* bases = GetSequence((hashes_collection.Get(id))->refBegins[0], special);
    //al->reference.assign(bases, (hashes_collection.Get(id))->length);
    //al->query.assign(bases, (hashes_collection.Get(id))->length);
    

    return true;
  }
}

void Aligner::GetTargetRefRegion(const int& read_length, const int& hash_begin, 
    const bool& special, int* begin, int* end) {
  uint32_t pos = hash_begin;
  int32_t ref_id = 0;
  if (special)
    SR_GetRefFromSpecialPos(special_ref_view_, &ref_id, &pos, reference_header_, reference_special_, hash_begin);

  int seq_length = special ? special_ref_view_->seqLen : reference_->seqLen;

  int forward_shift;
  if (pos < static_cast<uint32_t>(read_length)) forward_shift = pos;
  else forward_shift = read_length;

  int backward_shift;
  if ((pos + read_length + 1) > seq_length) backward_shift = seq_length - pos - 1;
  else backward_shift = read_length;
  
  *begin = hash_begin - forward_shift;
  *end   = hash_begin + backward_shift;
}

inline const char* Aligner::GetSequence(const size_t& start, const bool& special) const {
  if (special)
    return (reference_special_->sequence + start);
  else
    return (reference_->sequence + start);
}
