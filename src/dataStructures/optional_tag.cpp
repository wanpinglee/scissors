#include <list>
#include <string>
#include <iostream>
#include <outsources/samtools/bam.h>

using std::string;

namespace {
inline void ExtenedBamDataBuffer(
    const int& needed_size, 
    bam1_t* bam_record) {
  bam_record->m_data = needed_size;
  kroundup32(bam_record->m_data);
  uint8_t* ptr = bam_record->data;
  bam_record->data = (uint8_t*)calloc(bam_record->m_data, sizeof(uint8_t));
  memcpy(bam_record->data, ptr, bam_record->data_len);
  free(ptr);
}
}

namespace OptionalTag {
void AddOptionalTags(const bam1_t& anchor, std::list<bam1_t*> alignments) {
  const uint8_t* rg_ptr = bam_aux_get(&anchor, "RG");
  const char* rg = bam_aux2Z(rg_ptr);
  int rg_len = strlen(rg);

  for (std::list<bam1_t*>::iterator ite = alignments.begin();
       ite != alignments.end(); ++ite) {
    int ori_len = (*ite)->data_len;
    (*ite)->data_len += (3 + rg_len + 1);
    (*ite)->l_aux += (3 + rg_len + 1);
    if ((*ite)->m_data < (*ite)->data_len)
      ExtenedBamDataBuffer((*ite)->data_len, *ite);

    (*ite)->data[ori_len] = 'R';
    (*ite)->data[ori_len + 1] = 'G';
    (*ite)->data[ori_len + 2] = 'Z';
    memcpy((*ite)->data + ori_len + 3, rg, rg_len);
    (*ite)->data[ori_len + 3 + rg_len] = '\0';
  }
}

} // namespace OptionalTags
