#include <list>
#include <string>
#include <outsources/samtools/bam.h>

using std::string;

namespace OptionalTag {
void AddOptionalTags(const bam1_t& anchor, std::list<bam1_t*> alignments) {
  const uint8_t* rg_ptr = bam_aux_get(&anchor, "RG");
  const char* rg = bam_aux2Z(rg_ptr);
  int rg_len = strlen(rg);

  for (std::list<bam1_t*>::iterator ite = alignments.begin();
       ite != alignments.end(); ++ite) {
    bam_aux_append(*ite, "RG", 'Z', rg_len, (uint8_t*) rg);
  }
}

} // namespace OptionalTags
