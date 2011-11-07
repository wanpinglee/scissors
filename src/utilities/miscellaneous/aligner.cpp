#include "aligner.h"

#include "utilities/bam/bam_utilities.h"
#include "utilities/miscellaneous/hashes_collection.h"

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

    SR_SetStrand(query_region->pOrphan, SR_FORWARD);
    query_region->isOrphanInversed = FALSE;
  } // end if-else
  
}

Aligner::Aligner(const SR_Reference* reference, 
                 const SR_InHashTable* hash_table,
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


void Aligner::AlignCandidate(SR_BamInStreamIter* al_ite,
                             vector<bam1_t*>* alignments) {

    while (SR_QueryRegionLoadPair(query_region_, al_ite) == SR_OK) {
      
      if (query_region_->algnType != SR_UNIQUE_ORPHAN) continue;
      
      const bool is_anchor_forward = !bam1_strand(query_region_->pAnchor);
      //printf("%c\n", is_anchor_forward ? 'T' : 'F');
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
        
	//printf("%c\t%c\t%c\n",region_type.upstream ? 'T' : 'F', region_type.sequence_inverse ? 'T' : 'F', region_type.sequence_complement ? 'T' : 'F');
	// Reverse or complement the sequence if necesary
        SetTargetSequence(region_type, query_region_);
        HashRegionTableInit(hashes_, read_length);
        SR_QueryRegionSetRange(query_region_, &hash_length_, reference_->seqLen,
                               region_type.upstream ? SR_DOWNSTREAM : SR_UPSTREAM);
	HashRegionTableLoad(hashes_, hash_table_, query_region_);
	hashes_collection.Init(*(hashes_->pBestCloseRegions));
	//hashes_collection.Print();
	//char a; scanf("%c",&a);
	hashes_collection.SortByLength();
	//hashes_collection.Print();
	//scanf("%c",&a);

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

	  //printf("%s\n%s\n", al.reference.c_str(), al.query.c_str());
	  //alignments->push_back(al);
	  bam1_t* al_bam;
	  al_bam = bam_init1(); // Thread.cpp will free it
	  BamUtilities::ConvertAlignmentToBam1(al, *query_region_->pOrphan, al_bam);
	  alignments->push_back(al_bam);
	}

	//bam_write1( files.bam_writer, vars.query_region->pAnchor );
	//bam_write1( files.bam_writer, vars.query_region->pOrphan );
	//
      } // end while
    } // end while

    al_ite = NULL;
}
