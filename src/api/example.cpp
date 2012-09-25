#include <string>
#include <iostream>

#include "api.h"

using std::string;
using std::cout;
using std::endl;

const string kBamfile = "SARS.bam";
const string kFasta   = "SARS.reference.fa";

bool Prepare(string* fasta_seq, bam1_t* anchor, bam1_t* target, bam_header_t* header);
bool PrintAlignment(const bam1_t& al, const bam_header_t& header);

int main() {
  string fasta_seq; // reference bases are stored here
  bam1_t *anchor, *target; // anchor: mapped mate; target: unmapped mate
  anchor = bam_init1();
  target = bam_init1();
  bam_header_t* header = bam_header_init();

  // ====================================
  // Load reference bases from fasta file 
  // and anchor and target from bam file
  // ====================================
  Prepare(&fasta_seq, anchor, target, header);

  // ====================
  // Hash reference bases
  // ====================
  ReferenceHasher ref_hasher(fasta_seq.c_str());
  ref_hasher.Load();
  const SR_Reference* ref_ref    = ref_hasher.GetReference();
  const SR_InHashTable* ref_hash = ref_hasher.GetHashTable();
  // [NOTICE]: If you like to detect special insertions
  // SpecialHasher sp_hasher(special_fasta_filename);
  // sp_hasher.Load();
  // const SR_Reference* sp_ref = sp_hasher.GetReference();
  // const SR_RefHeader* sp_header = sp_hasher.GetReferenceHeader();
  // const SR_InHashTable* sp_hash = sp_hasher.GetHashTable();

  // ============================
  // Declare a split-read aligner
  // ============================
  Scissors::Aligner aligner;
  aligner.SetReference(ref_ref, ref_hash);
  // [NOTICE]: If you like to detect special insertions
  // Scissors::Technology tech;
  // tech = Scissors::TECH_ILLUMINA;
  // aligner.SetReference(ref_ref, ref_hash, tech, sp_ref, sp_hash, sp_header);

  // ===========================
  // Split-read align the target
  // ===========================
  Scissors::TargetEvent event; // default TargetEvent
  Scissors::TargetRegion region; // default TargetRegion
  region.fragment_length = 1000; // set the fragment length
  Scissors::AlignmentFilter filter; // default AlignmentFilter
  vector<bam1_t*> result; // resultant alignments are sotred here
                          // [NOTICE]: bamt_1* has to be freed
  // [NOTICE]: If you like to detect special insertions
  // event.special_insertion = true;
  // *** Aligning ***
  aligner.AlignCandidate(event, region, filter, *anchor, *target, &result);

  // ==============================
  // Print the resultant alignments
  // ==============================
  cout << "[RESULTS]" << endl;
  for (unsigned int i = 0; i < result.size(); ++i)
    PrintAlignment(*result[i], *header);
  
  // ===============
  // Clean up memory
  // ===============
  bam_destroy1(anchor);
  bam_destroy1(target);
  bam_header_destroy(header);

  for (unsigned int i = 0; i < result.size(); ++i)
    bam_destroy1(result[i]);


  return 0;
}

// =========
// functions
// =========
bool Prepare(string* fasta_seq, bam1_t* anchor, bam1_t* target, bam_header_t* header) {
  // Load references fasta
  FastaReference ref_reader;
  ref_reader.open(kFasta.c_str());
  string ref_name;
  if (ref_reader.index->sequenceNames.size() == 0) return false;
  else ref_name = *(ref_reader.index->sequenceNames.begin());
  *fasta_seq = ref_reader.getSequence(ref_name);

  // Load bam records
  bamFile bam_reader = bam_open(kBamfile.c_str(), "r");
  header = bam_header_read(bam_reader);
  bam_read1(bam_reader, target); // unmapped mate
  bam_read1(bam_reader, anchor); // mapped mate
  //anchor->core.tid = 0;
  //anchor->core.pos = 23539;
  //anchor->core.flag = 147;
  bam_close(bam_reader);

  return true;
}

bool PrintAlignment(const bam1_t& al, const bam_header_t& header) {
  cout << bam1_qname(&al) << "\t"
       << al.core.flag << "\t"
       << "chrID:" << al.core.tid << "\t"
       << al.core.pos << "\t"
       << al.core.qual << "\t";
  
  uint32_t* cigar = bam1_cigar(&al);
  for (unsigned int i = 0; i < al.core.n_cigar; ++i) {
    cout << (cigar[i] >> 4);
    int op = cigar[i] & 0x0000000f;
    switch (op) {
      case 0: cout << 'M'; break;
      case 1: cout << 'I'; break;
      case 2: cout << "D"; break;
      case 3: cout << 'N'; break;
      case 4: cout << 'S'; break;
      case 5: cout << 'H'; break;
      case 6: cout << 'P'; break;
      case 7: cout << '='; break;
      case 8: cout << 'X'; break;
      default: cout << '?';
    }
  }
  cout << "\t";

  cout << al.core.mtid << "\t"
       << al.core.mpos << "\t";

  for (int i = 0; i < al.core.l_qseq; ++i) {
    int base = bam1_seqi(bam1_seq(&al), i);
    switch (base) {
      case 1: cout << 'A'; break;
      case 2: cout << 'C'; break;
      case 4: cout << 'G'; break;
      case 8: cout << 'T'; break;
      case 15: cout << 'N'; break;
      default: cout << '?';
    }
  }
  cout << "\t";

  char* ptr = (char*) bam1_qual(&al);
  for (int i = 0; i < al.core.l_qseq; ++i) {
    cout << (char) ((*ptr) + 33);
    ++ptr;
  }
  cout << "\t";
       
  cout << bam1_aux(&al) << endl;
  return true;
}
