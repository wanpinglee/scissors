#ifndef _ALIGNER_H_
#define _ALIGNER_H_

extern "C" {
#include "dataStructures/SR_QueryRegion.h"
#include "outsources/samtools/bam.h"
#include "utilities/bam/SR_BamInStream.h"
#include "utilities/hashTable/SR_HashRegionTable.h"
#include "utilities/hashTable/SR_InHashTable.h"
#include "utilities/hashTable/SR_Reference.h"
}

#include "dataStructures/anchor_region.h"
#include "dataStructures/search_region_type.h"
#include "dataStructures/technology.h"
#include "utilities/smithwaterman/ssw_cpp.h"

namespace Scissors {

struct Alignment;
struct AlignmentFilter;
struct TargetEvent;
struct TargetRegion;
class HashesCollection;

class Aligner {
 public:
  Aligner();
  Aligner(const SR_Reference*   reference, 
          const SR_InHashTable* hash_table,
	  const SR_Reference*   reference_special,
	  const SR_InHashTable* hash_table_special,
	  const SR_RefHeader*   reference_header,
	  const Technology&     technology);
  ~Aligner();
 
  // @function:
  //     Set reference and hash table
  // @return:
  //     always true
  //
  // @params:
  //     If not detect special references, then reference_special,
  //     hash_table_special, and reference_header are not necessary to set.
  bool SetReference(const SR_Reference*   reference,
                    const SR_InHashTable* hash_table,
		    const Technology&     technology         = TECH_ILLUMINA,
		    const SR_Reference*   reference_special  = NULL,
		    const SR_InHashTable* hash_table_special = NULL,
		    const SR_RefHeader*   reference_header   = NULL);

  // [NOTICE] Users may not use this function.
  // @function:
  //     Aligns the orphan bam1_t in query_region_
  //
  // @return:
  //     false: when (SR_Reference*) reference is not set 
  //            or (Technology) technology is TECH_NONE;
  //     true: otherwise
  //
  // @params:
  //     target_event------the events that are gonna be checked.
  //     target_region-----*fragment_length(, local_window_size,
  //                       and discovery_window_size are set here.
  //     alignment_filter--the filter for split-read alignment.
  //     al_ite------------all target pairs are stored here;
  //                       it'll be set to NULL before exiting the function.
  //     alignments--------all obtained split-read alignments are stored here;
  //                       [NOTICE] Users should free bam1_t in the vector.
  bool AlignCandidate(const TargetEvent&     target_event,
                      const TargetRegion&    target_region,
                      const AlignmentFilter& alignment_filter,
                      SR_BamInStreamIter*    al_ite, 
                      vector<bam1_t*>*       alignments);

  // @function:
  //     Aligns the orphan bam1_t in query_region_
  //
  // @return:
  //     false: when (SR_Reference*) reference is not set
  //            or (Technology) technology is TECH_NONE;
  //     true: otherwise
  //
  // @params:
  //     target_event------the events that are gonna be checked
  //     target_region-----*fragment_length(, local_window_size, 
  //                       and discovery_window_size are set here.
  //     alignment_filter--the filter for split-read alignment.
  //     anchor------------anchor (mapped) mate.
  //                       [NOTICE] Only tid, pos, and flag are used.
  //     target------------target (unmapped) mate which is the split-read target;
  //                       [NOTICE] tid, flag, bases, and qualities are needed.
  //     alignments--------all obtained split-read alignments are stored here;
  //                       [NOTICE] Users should free bam1_t in the vector.
  bool AlignCandidate(const TargetEvent&     target_event,
                      const TargetRegion&    target_region,
                      const AlignmentFilter& alignment_filter,
		      const bam1_t&          anchor,
		      const bam1_t&          target,
		      vector<bam1_t*>*       alignments);
  // @function:
  //     Change the Smith-Waterman score matrix 
  //     for the medium-sized-indel ssw aligner
  bool ResetIndelSmithWatermanScore(const uint8_t& match_score      = 10,
                                    const uint8_t& mismatch_penalty      = 20,
	                            const uint8_t& gap_opening_penalty   = 20,
		                    const uint8_t& gap_extending_penalty = 1);

  // @function:
  //     Change the Smith-Waterman score matrix for the local ssw aligner
  bool ResetLocalSmithWatermanScore(const uint8_t& match_score      = 2,
                                    const uint8_t& mismatch_penalty      = 2,
	                            const uint8_t& gap_opening_penalty   = 3,
		                    const uint8_t& gap_extending_penalty = 1);
 private:
  SearchRegionType search_region_type_;
  AnchorRegion     anchor_region_;
  Technology technology_;

  const SR_Reference*   reference_;
  const SR_InHashTable* hash_table_;
  const SR_Reference*   reference_special_;
  const SR_InHashTable* hash_table_special_;
  const SR_RefHeader*   reference_header_;
  SR_QueryRegion*       query_region_;
  HashRegionTable*      hashes_;
  HashRegionTable*      hashes_special_;
  SR_SearchArgs         hash_length_;
  SR_RefView*           special_ref_view_;

  StripedSmithWaterman::Aligner stripe_sw_indel_;
  StripedSmithWaterman::Aligner stripe_sw_normal_;

  void LoadRegionType(const bam1_t& anchor);
  const char* GetSequence(const size_t& start, const bool& special) const;
  bool GetAlignment(const HashesCollection& hashes_collection, 
                    const unsigned int& id, const bool& special, const int& read_length,
                    const char* read_seq, StripedSmithWaterman::Alignment* al);
  void GetTargetRefRegion(const int& extend_length, const int& hash_begin,
                          const bool& special, int* begin, int* end);
  void Align(const TargetEvent& target_event,
             const TargetRegion& target_region,
             const AlignmentFilter& alignment_filter,
	     const SR_QueryRegion* query_region_,
	     vector<bam1_t*>* alignments);
  bool SearchLocalPartial(const TargetRegion& target_region,
                          const AlignmentFilter& alignment_filter,
			  StripedSmithWaterman::Alignment* local_al);
  bool SearchSpecialReference(const TargetRegion& target_region,
                              const AlignmentFilter& alignment_filter,
			      StripedSmithWaterman::Alignment* special_al);
  bool SearchMediumIndel(const TargetRegion& target_region,
                         const AlignmentFilter& alignment_filter,
                         StripedSmithWaterman::Alignment* ssw_al);

  Aligner (const Aligner&);
  Aligner& operator= (const Aligner&);
};
} // namespace
#endif // _ALIGNER_H_
