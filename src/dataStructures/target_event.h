#ifndef DATASTRUCTURES_TARGET_EVENT_H_
#define DATASTRUCTURES_TARGET_EVENT_H_

namespace Scissors {
struct TargetEvent {
  int  local_region;
  bool special_insertion;
  bool medium_sized_indel;

  TargetEvent()
      : special_insertion(false)
      , medium_sized_indel(true)
  {}
};
} // namespace
#endif // DATASTRUCTURES_TARGET_EVENT_H_
