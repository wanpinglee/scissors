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
#include "utilities/smithwaterman/BandedSmithWaterman.h"
#include "utilities/smithwaterman/ssw_cpp.h"

namespace Scissors {

struct Alignment;
struct AlignmentFilter;
struct TargetEvent;
struct TargetRegion;
class HashesCollection;

class Aligner {
 public:
  Aligner(const SR_Reference*   reference, 
          const SR_InHashTable* hash_table,
	  const SR_Reference*   reference_special,
	  const SR_InHashTable* hash_table_special,
	  const SR_RefHeader*   reference_header,
	  const int&            fragment_length);
  ~Aligner();
  // @function Aligns the orphan bam1_t in query_region_
  //           [NOTICE] Users may not use this function.
  // @param  detect_special    detects second partial in special references
  // @param  alignment_filter  the filter for split-read alignment
  // @param  al_ite            all target pairs are stored here
  // @param  alignments        all obtained split-read alignments are stored here
  //         [NOTICE]          Users should free bam1_t in the vector
  void AlignCandidate(const TargetEvent&     target_event,
                      const TargetRegion&    target_region,
                      const AlignmentFilter& alignment_filter,
                      SR_BamInStreamIter*    al_ite, 
                      vector<bam1_t*>*       alignments);

  // @function Aligns the orphan bam1_t in query_region_
  // @param  detect_special    detects second partial in special references
  // @param  alignment_filter  the filter for split-read alignment
  // @param  anchor            anchor (mapped) mate
  // @param  target            target (unmapped) mate which is the split-read target
  // @param  alignments        all obtained split-read alignments are stored here
  //         [NOTICE]          Users should free bam1_t in the vector
  void AlignCandidate(const TargetEvent&     target_event,
                      const TargetRegion&    target_region,
                      const AlignmentFilter& alignment_filter,
		      const bam1_t&          anchor,
		      const bam1_t&          target,
		      vector<bam1_t*>*       alignments);
 private:
  SearchRegionType search_region_type_;
  AnchorRegion     anchor_region_;

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

  CBandedSmithWaterman  banded_sw_aligner_;
  StripedSmithWaterman::Aligner stripe_sw_aligner_;

  void LoadRegionType(const bam1_t& anchor);
  const char* GetSequence(const size_t& start, const bool& special) const;
  bool GetAlignment(const HashesCollection& hashes_collection, 
                    const unsigned int& id, const bool& special, const int& read_length,
                    const char* read_seq, Alignment* al);
  void GetTargetRefRegion(const int& extend_length, const int& hash_begin,
                          const bool& special, int* begin, int* end);
  void Align(const TargetEvent& target_event,
             const TargetRegion& target_region,
             const AlignmentFilter& alignment_filter,
	     const SR_QueryRegion* query_region_,
	     vector<bam1_t*>* alignments);
  bool SearchSpecialReference(const AlignmentFilter& alignment_filter,
                              Alignment* local_al, Alignment* special_al);
  void SearchLocalRegion(const TargetRegion& target_region, 
                         StripedSmithWaterman::Alignment* ssw_al);

  Aligner (const Aligner&);
  Aligner& operator= (const Aligner&);
};
} // namespace
#endif // _ALIGNER_H_
