#ifndef UTILITIES_MISCELLANEOUS_ALIGNMENT_COLLECTION_H_
#define UTILITIES_MISCELLANEOUS_ALIGNMENT_COLLECTION_H_

#include <vector>

#include "utilities/smithwaterman/ssw_cpp.h"

namespace Scissors {

struct Alignment;
struct TargetEvent;

// Collects split-read alignments that support different events
class AlignmentCollection {
 public:
  AlignmentCollection();
  ~AlignmentCollection();

  // PushANewEvent:
  // @function:
  //     A new event is supported by split-read alignment and
  //     the function will push the event in the container.
  //     The alignments supporting the event should be pushed by
  //     PushAlignment function.
  // @params:
  //     event--the type of the event
  void PushANewEvent(const TargetEvent& event);

  // PushAlignment:
  // @function:
  //     Push an SSW alignment that supports the pushed event
  //     [NOTICE] The alignment corresponds to the latest event
  void PushAlignment(const StripedSmithWaterman::Alignment& al);

  // PushAlignment:
  // @function:
  //     Push a common alignment that supports the pushed event
  //     [NOTICE] The alignment corresponds to the latest event
  void PushAlignment(const Alignment& al);

  // GetMostConfidentEvent:
  // @function:
  //     Seek the most confident event which means the event has 
  //     the best alignment support, and then report event type
  //     and its corresponding alignments.
  // @params:
  //     event-------the most confident event
  //     ssw_als-----the ssw alignments supporting the most confident event
  //     common_als--the common alignments supporting the most confident event
  void GetMostConfidentEvent(
      TargetEvent* event, 
      std::vector <StripedSmithWaterman::Alignment*>* ssw_als,
      std::vector <Alignment*>* common_als);

  // Clear:
  // @function:
  //     Clear the container
  void Clear();

  struct EventAlignment{
    TargetEvent* event;
    std::vector <StripedSmithWaterman::Alignment*> ssw_als;
    std::vector <Alignment*> common_als;
  };
 
 private:

  std::vector <EventAlignment> container_;

  AlignmentCollection (const AlignmentCollection&);
  AlignmentCollection& operator= (const AlignmentCollection&);

}; // class AlignmentCollection
} // namespace Scissors

#endif // UTILITIES_MISCELLANEOUS_ALIGNMENT_COLLECTION_H_
