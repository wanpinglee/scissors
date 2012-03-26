#ifndef DATASTRUCTURES_OPTIONAL_TAG_H_
#define DATASTRUCTURES_OPTIONAL_TAG_H_

#include <list>
#include <outsources/samtools/bam.h>

namespace OptionalTag {

bool AddOptionalTags (const bam1_t& anchor, std::list<bam1_t*>alignments);

} // namespace OptionalTags

#endif // DATASTRUCTURES_OPTIONAL_TAG_H_
