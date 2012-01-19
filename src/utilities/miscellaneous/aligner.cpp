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
    , hash_table_(hash_table) {
  query_region_ = SR_QueryRegionAlloc();
  hashes_       = HashRegionTableAlloc();

  hash_length_.fragLen = fragment_length;
  hash_length_.closeRange = 2000;
  hash_length_.farRange   = 10000;
}

Aligner::~Aligner() {
  SR_QueryRegionFree(query_region_);
  HashRegionTableFree(hashes_);
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
      
      // store the anchor in output bam
      bam1_t* al_bam_anchor;
      al_bam_anchor = bam_init1(); // Thread.cpp will free it
      al_bam_anchor = bam_copy1(al_bam_anchor, query_region_->pAnchor);
      alignments->push_back(al_bam_anchor);

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
	/*
	int hashes_count = hashes_collection.GetSize();
	if (hashes_count == 0) continue;
	if ((hashes_collection.Get(hashes_count-1))->numPos != 0) {
	  Alignment al;
	  al.reference_begin = (hashes_collection.Get(hashes_count-1))->refBegins[0];
	  al.reference_end   = (hashes_collection.Get(hashes_count-1))->refBegins[0] + (hashes_collection.Get(hashes_count-1))->length - 1;
	  al.query_begin     = (hashes_collection.Get(hashes_count-1))->queryBegin;
	  al.query_end       = (hashes_collection.Get(hashes_count-1))->queryBegin + (hashes_collection.Get(hashes_count-1))->length - 1;
	  al.is_seq_inverse  = region_type.sequence_inverse;
	  al.is_seq_complement = region_type.sequence_complement;
	  const char* bases = GetSequence((hashes_collection.Get(hashes_count-1))->refBegins[0]);
	  al.reference.assign(bases, (hashes_collection.Get(hashes_count-1))->length);
	  al.query.assign(bases, (hashes_collection.Get(hashes_count-1))->length);
	  
	  bam1_t* al_bam;
	  al_bam = bam_init1(); // Thread.cpp will free it
	  BamUtilities::ConvertAlignmentToBam1(al, *query_region_->pOrphan, al_bam);
	  alignments->push_back(al_bam);
	}
	*/
      } // end while
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
