#ifndef DATASTRUCTURES_OPTIONAL_TAG_H_
#define DATASTRUCTURES_OPTIONAL_TAG_H_

#include <list>
#include <outsources/samtools/bam.h>

namespace Scissors {
namespace OptionalTag {

bool AddOptionalTags (const bam1_t& anchor, bam1_t* alignments);
bool AddOptionalTags (const bam1_t& anchor, std::list<bam1_t*>alignments);

} // namespace OptionalTags
} // namespace Scissors

#endif // DATASTRUCTURES_OPTIONAL_TAG_H_
