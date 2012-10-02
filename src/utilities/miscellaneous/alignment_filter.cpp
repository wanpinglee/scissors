#include "alignment_filter.h"

#include <iostream>
#include "dataStructures/alignment.h"
#include "utilities/smithwaterman/ssw_cpp.h"

namespace Scissors {
namespace AlignmentFilterApplication {
void TrimAlignment(
    const AlignmentFilter& filter, 
    StripedSmithWaterman::Alignment* al) {
  // find the sublist whose scores is max
  int maxi = 0, max = 0, start = 0, length = 0, best_start = 0, best_length = 0;
  for (unsigned int i = 0; i < al->cigar.size(); ++i) {
    uint32_t op = al->cigar[i] & 0x0000000f;
    uint32_t op_length = al->cigar[i] >> 4;
    // get score
    int score = 0;
    switch(op) {
      case 0: // M
        score = filter.trimming_match_score * op_length;
	break;
      case 1: // I
      case 2: // D
        score = filter.trimming_gap_score * op_length;
	break;
      default:
        break;
    } // end of switch
    
    // find max
    if ((maxi + score) > 0) {
      maxi += score;
      ++length;
    } else {
      maxi = 0;
      start = i + 1;
      length = 0;
    }
    if (maxi > max) {
      max = maxi;
      best_start  = start;
      best_length = length;
    }

    //std::cerr << maxi << " " << start << " " << length << std::endl;
  } // end of for

  //std::cerr << max << " " << best_start << " " << best_length << std::endl;

  // trim alignment
  if (best_length == 0) {
    al->cigar.clear();
    al->cigar_string.clear();
  } else if (static_cast<unsigned int>(best_length) == al->cigar.size()){
    // nothing
  } else {
    int read_clip1  = 0;
    int ref_clip1   = 0;
    int query_clip1 = 0;
    int mismatch1   = 0;
    // Calculate the beginning soft clip
    for (unsigned int i = 0; i < static_cast<unsigned int>(best_start); ++i) {
      uint32_t op = al->cigar[i] & 0x0000000f;
      uint32_t op_length = al->cigar[i] >> 4;
      switch(op) {
        case 0: // M
	  read_clip1  += op_length;
	  ref_clip1   += op_length;
	  query_clip1 += op_length;
	  break;
	case 4: // S
	  read_clip1 += op_length;
	  break;
	case 1: // I
	  read_clip1  += op_length;
	  query_clip1 += op_length;
	  mismatch1   += op_length;
	  break;
	case 2: // D
	  mismatch1 += op_length;
	  ref_clip1 += op_length;
	  break;
	default:
	  break;
      }
    } // end of for

    //std::cerr << "begin clip: " << read_clip1 << std::endl;

    // Calculate the ending clip
    int read_clip2  = 0;
    int ref_clip2   = 0;
    int query_clip2 = 0;
    int mismatch2   = 0;
    for (unsigned int i = best_start + best_length; i < al->cigar.size(); ++i) {
      uint32_t op = al->cigar[i] & 0x0000000f;
      uint32_t op_length = al->cigar[i] >> 4;
      switch(op) {
        case 0: // M
	  read_clip2  += op_length;
	  ref_clip2   += op_length;
	  query_clip2 += op_length;
	  break;
	case 4: // S
	  read_clip2 += op_length;
	  break;
	case 1: // I
	  read_clip2  += op_length;
	  query_clip2 += op_length;
	  mismatch2   += op_length;
	  break;
	case 2: // D
	  mismatch2 += op_length;
	  ref_clip2 += op_length;
	  break;
	default:
	  break;
      }
    } // end of for

    //std::cerr << "end clip: " << read_clip2 << std::endl;
    
    // Modify cigar vector
    uint32_t new_cigar = 0;
    if (read_clip2 > 0) {
      al->cigar.erase(al->cigar.begin() + best_start + best_length, al->cigar.end());
      new_cigar = (read_clip2 << 4 ) | 0x04;
      al->cigar.push_back(new_cigar);
      al->ref_end    -= ref_clip2;
      al->query_end  -= query_clip2;
      al->mismatches -= mismatch2;
    }
    new_cigar = 0;
    if (read_clip1 > 0) {
      al->cigar.erase(al->cigar.begin(), al->cigar.begin() + best_start);
      new_cigar = (read_clip1 << 4 ) | 0x04;
      al->cigar.insert(al->cigar.begin(), new_cigar);
      al->ref_begin   += ref_clip1;
      al->query_begin += query_clip1;
      al->mismatches  -= mismatch1;
    }

    if (al->mismatches < 0) // it should not happen
      al->mismatches = 0;
  }
}

void TrimAlignment(
    const AlignmentFilter& filter,
    Alignment* al){

  // find the sublist whose scores is max
  int maxi = 0, max = 0, start = 0, length = 0, best_start = 0, best_length = 0;
  for (unsigned int i = 0; i < al->reference.size(); ++i) {
    // get score
    int score = 0;
    if ((al->reference[i] == AlignmentConstant::GAP) || (al->query[i] == AlignmentConstant::GAP)) {
      score = filter.trimming_gap_score;
    } else {
      if (al->reference[i] != al->query[i]) score = filter.trimming_mismatch_score;
      else score = filter.trimming_match_score;
    }

    // find max
    if ((maxi + score) > 0) {
      maxi += score;
      ++length;
    } else {
      maxi = 0;
      start = i + 1;
      length = 0;
    }
    if (maxi > max) {
      max = maxi;
      best_start = start;
      best_length = length;
    }
  }

  // trim alignment
  if ((best_length > 0) 
  && (static_cast<unsigned int>(best_length) < al->reference.size())) {
    int end = best_start + best_length - 1;
    //sanity check
    if (end > (al->reference.size() - 1)) {
      //nothing
    } else {
      // trim from the tail
      if (end < (al->reference.size() - 1)) {
        for (unsigned int i = end + 1; i < al->reference.size(); ++i) {
	  if (al->reference[i] != AlignmentConstant::GAP) --(al->reference_end);
	  if (al->query[i] != AlignmentConstant::GAP) --(al->query_end);
	  if (al->reference[i] != al->query[i]) --(al->num_mismatches);
	}
	int size = al->reference.size();
	al->reference.erase(end + 1, size - end - 1);
	al->query.erase(end + 1, size - end - 1);
      }
      
      // trim from the begin
      if (best_start > 0) {
        for (int i = 0; i < best_start; ++i) {
	  if (al->reference[i] != AlignmentConstant::GAP) ++(al->reference_begin);
	  if (al->query[i] != AlignmentConstant::GAP) ++(al->query_begin);
	  if (al->reference[i] != al->query[i]) --(al->num_mismatches);
	}
        al->reference.erase(0, best_start);
        al->query.erase(0, best_start);
      }
    }
  } else {
    // nothing
  }
} // TrimAlignment

} //namespace AlignmentFilter
} //namespaceScissors
