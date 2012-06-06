#include "utilities/miscellaneous/alignment_collection.h"

#include <algorithm>

#include "dataStructures/target_event.h"
#include "dataStructures/alignment.h"

using std::vector;

namespace Scissors {
namespace {
bool SortEventAlignment(const AlignmentCollection::EventAlignment& ea1, const AlignmentCollection::EventAlignment& ea2) {
  int ea1_base = 0;
  for (vector <StripedSmithWaterman::Alignment*>::const_iterator ite = ea1.ssw_als.begin();
      ite != ea1.ssw_als.end(); ++ite) {
    ea1_base += (*ite)->query_end - (*ite)->query_begin;
  }
  for (vector <Alignment*>::const_iterator ite = ea1.common_als.begin();
      ite != ea1.common_als.end(); ++ite) {
    ea1_base += (*ite)->query_end - (*ite)->query_begin;
  }

  int ea2_base = 0;
  for (vector <StripedSmithWaterman::Alignment*>::const_iterator ite = ea2.ssw_als.begin();
      ite != ea2.ssw_als.end(); ++ite) {
    ea2_base += (*ite)->query_end - (*ite)->query_begin;
  }
  for (vector <Alignment*>::const_iterator ite = ea2.common_als.begin();
      ite != ea2.common_als.end(); ++ite) {
    ea2_base += (*ite)->query_end - (*ite)->query_begin;
  }

  return ea1_base < ea2_base;

}
} // unnamed namespace

AlignmentCollection::AlignmentCollection() {
}

AlignmentCollection::~AlignmentCollection() {
  Clear();
}

// PushANewEvent:
// @function:
//     A new event is supported by split-read alignment and
//     the function will push the event in the container.
//     The alignments supporting the event should be pushed by
//     PushAlignment function.
// @params:
//     event--the type of the event
void AlignmentCollection::PushANewEvent(const TargetEvent& event) {
  EventAlignment new_result;
  new_result.event = event;
  container_.push_back(new_result);
}

void AlignmentCollection::PushAlignment(
    const StripedSmithWaterman::Alignment& al) {
  if (container_.empty()) {
    return;
  } else {
    // Gets the lastest event
    vector<EventAlignment>::reverse_iterator ite = container_.rbegin();
    ite->ssw_als.push_back((StripedSmithWaterman::Alignment*)&al);
  }
}

void AlignmentCollection::PushAlignment(const Alignment& al) {
  if (container_.empty()) {
    return;
  } else {
    // Gets the lastest event
    vector<EventAlignment>::reverse_iterator ite = container_.rbegin();
    ite->common_als.push_back((Alignment*)&al);
  }
}

void AlignmentCollection::GetMostConfidentEvent(
    TargetEvent* event,
    std::vector <StripedSmithWaterman::Alignment*>* ssw_als,
    std::vector <Alignment*>* common_als) {

  if (container_.empty()) {
    event->special_insertion  = false;
    event->medium_sized_indel = false;
    ssw_als->clear();
    common_als->clear();
  } else {
    std::sort(container_.begin(), container_.end(), SortEventAlignment);
    vector<EventAlignment>::reverse_iterator ite = container_.rbegin();
    *event      = ite->event;
    *ssw_als    = ite->ssw_als;
    *common_als = ite->common_als;
  }
}

// @function:
//     Clear the container
void AlignmentCollection::Clear() {
  for (vector<EventAlignment>::iterator ite = container_.begin(); 
      ite != container_.end(); ++ite) {
    ite->ssw_als.clear();
    ite->common_als.clear();
  }

  container_.clear();
}
} // namespace Scissors
