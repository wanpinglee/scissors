#include "aligner.h"
#include "utilities/bam/bam_utilities.h"
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
}

Aligner::Aligner(const SR_Reference* reference, 
                 const SR_InHashTable* hash_table,
		 const SR_Reference* reference_special,
		 const SR_InHashTable* hash_table_special,
		 const int& fragment_length)
    : reference_(reference)
    , hash_table_(hash_table)
    , reference_special_(reference_special)
    , hash_table_special_(hash_table_special) {
  query_region_   = SR_QueryRegionAlloc();
  hashes_         = HashRegionTableAlloc();
  hashes_special_ = HashRegionTableAlloc();

  hash_length_.fragLen = fragment_length;
  hash_length_.closeRange = 2000;
  hash_length_.farRange   = 10000;
}

Aligner::~Aligner() {
  SR_QueryRegionFree(query_region_);
  HashRegionTableFree(hashes_);
  HashRegionTableFree(hashes_special_);
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
                             SR_BamInStreamIter* al_ite,
                             vector<bam1_t*>* alignments) {

    while (SR_QueryRegionLoadPair(query_region_, al_ite) == SR_OK) {
      // TODO@WP: it may be removed later
      if (query_region_->algnType != SR_UNIQUE_ORPHAN) continue;
      
      const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
      // Convert 4-bit representive sequence into chars
      SR_QueryRegionLoadSeq(query_region_);
      int read_length = query_region_->pOrphan->core.l_qseq;

      LoadRegionType(*(query_region_->pAnchor));

      SearchRegionType::RegionType region_type;
      HashesCollection hashes_collection;
      HashesCollection hashes_collection_special;
      
      // store the anchor in output bam
      bam1_t* al_bam_anchor;
      al_bam_anchor = bam_init1(); // Thread.cpp will free it
      al_bam_anchor = bam_copy1(al_bam_anchor, query_region_->pAnchor);
      alignments->push_back(al_bam_anchor);

      // For MEI first
      search_region_type_.GetStandardType(is_anchor_forward, &region_type);
      SetTargetSequence(region_type, query_region_);
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
      hashes_collection_special.Print();

      // get best cover
      unsigned int best1, best2;
      bool best_pair_found = false;
      if (is_anchor_forward)
        best_pair_found = hashes_collection.GetBestCoverPair(&hashes_collection_special, &best1, &best2);
      else 
        best_pair_found = hashes_collection_special.GetBestCoverPair(&hashes_collection, &best2, &best1);
      Alignment al1, al2;
      if (best_pair_found) {
        GetAlignment(hashes_collection, best1, &al1);
	GetAlignment(hashes_collection_special, best2, &al2);
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
    Alignment* al) {
  if (static_cast<int>(id) >= hashes_collection.GetSize()) {
    return false;
  } else {
    al->reference_begin = (hashes_collection.Get(id))->refBegins[0];
    al->reference_end   = (hashes_collection.Get(id))->refBegins[0] + (hashes_collection.Get(id))->length - 1;
    al->query_begin     = (hashes_collection.Get(id))->queryBegin;
    al->query_end       = (hashes_collection.Get(id))->queryBegin + (hashes_collection.Get(id))->length - 1;
    const char* bases = GetSequence((hashes_collection.Get(id))->refBegins[0]);
    al->reference.assign(bases, (hashes_collection.Get(id))->length);
    al->query.assign(bases, (hashes_collection.Get(id))->length);

    return true;
  }
}
