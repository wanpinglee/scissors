#include "gtest/gtest.h"

#include <string>

extern "C" {
#include "outsources/samtools/bam.h"
#include "utilities/bam/SR_BamPairAux.h"
}

#include "utilities/hashTable/special_hasher.h"
#include "utilities/hashTable/reference_hasher.h"
#include "googletest/util/fasta_reader.h"
#include "utilities/miscellaneous/aligner.h"
#include "utilities/miscellaneous/alignment_filter.h"

using std::string;

const string path = "/share0/home/wanping/source/scissors/src/googletest/util/";

namespace {
void LoadSpecialReference(SpecialHasher* sp_hasher) {
  const string fasta = path + "moblist_19Feb2010.fa";
  sp_hasher->SetFastaName(fasta.c_str());
  ASSERT_TRUE(sp_hasher->Load());
}

void StoreAlignmentInBam(vector<bam1_t*>* alignments_bam,
                         bamFile* bam_writer) {
  for (unsigned int i = 0; i < alignments_bam->size(); ++i) {
    bam_write1(*bam_writer, (*alignments_bam)[i]);
  }
  for (unsigned int i = 0; i < alignments_bam->size(); ++i) {
    bam1_t* ptr = (*alignments_bam)[i];
    bam_destroy1(ptr);
  }

}

TEST(AlignerApi, Test) {
  // load special references; MEIs fasta
  SpecialHasher sp_hasher;
  // =====================================
  // build a hasher for special references
  // =====================================
  
  LoadSpecialReference(&sp_hasher);
  const SR_Reference* sp_ref        = sp_hasher.GetReference();
  const SR_RefHeader* sp_ref_header = sp_hasher.GetReferenceHeader();
  const SR_InHashTable* sp_hash     = sp_hasher.GetHashTable();

  ASSERT_TRUE(sp_ref != NULL);
  ASSERT_TRUE(sp_ref_header != NULL);
  ASSERT_TRUE(sp_hash != NULL);


  ASSERT_TRUE(sp_ref_header != NULL);
  ASSERT_TRUE(sp_hash != NULL);


  // load a chromosome; using chr20 for testing
  FastaReader fr;
  const string fasta = path + "hg19.fasta";
  fr.Open(fasta.c_str());
  string ref_name, ref_seq;
  // load until hitting chr20
  while (ref_name != "20") {
    if (!fr.LoadNextRead(&ref_name, &ref_seq)) break;
  }

  // ========================
  // build a hasher for chr20
  // ========================
  ReferenceHasher ref_hasher(ref_seq.c_str());
  ASSERT_TRUE(ref_hasher.Load());
  const SR_Reference* ref_ref    = ref_hasher.GetReference();
  const SR_InHashTable* ref_hash = ref_hasher.GetHashTable();

  ASSERT_TRUE(ref_ref != NULL);
  ASSERT_TRUE(ref_hash != NULL);
  fr.Close();
/*
  // ===============
  // bam preparation
  // ===============
  string bamoutput = path + "test.bam";
  string baminput = path + "NA18523.chr20.unampped.bam";
  bamFile bam_writer = bam_open(bamoutput.c_str(), "w");
  bamFile bam_reader = bam_open(baminput.c_str(), "r");
  bam_header_t* header = bam_header_init();
  header = bam_header_read(bam_reader);
  bam_header_write(bam_writer, header);
  bam_header_destroy(header);
  // end of preparation

  // ===================
  // declare the aligner
  // ===================
  Scissors::Technology tech = Scissors::TECH_ILLUMINA;
  Scissors::Aligner aligner(ref_ref, ref_hash, sp_ref, sp_hash, sp_ref_header, 200, tech);
  
  // =============================================
  // prepare alignment filter and output container
  // =============================================
  Scissors::AlignmentFilter alignment_filter;
  vector<bam1_t*> alignments_bam;

  bam1_t* anchor = bam_init1();
  bam1_t* target = bam_init1();
  int bam_status = 1;

  // =======================
  // Align and store results
  // =======================
  //while (bam_status > 0) {
  //  bam_status = bam_read1(bam_reader, anchor);
  //  bam_status = bam_read1(bam_reader, target);
  //  alignments_bam.clear();
  //  aligner.AlignCandidate(true, alignment_filter, *anchor, *target, &alignments_bam);
  //  StoreAlignmentInBam(&alignments_bam, &bam_writer);
  //}

  bam_destroy1(anchor);
  bam_destroy1(target);
  bam_close(bam_reader);
  bam_close(bam_writer);
*/
  // =================================
  // demo how to reuse ReferenceHasher
  // =================================

  ref_hasher.Clear();
  ref_hasher.Clear();
  ref_hasher.Clear();

  ref_hasher.SetSequence(ref_seq.c_str());
  ASSERT_TRUE(ref_hasher.Load());
  ref_hasher.Clear();
}

} //namespace
