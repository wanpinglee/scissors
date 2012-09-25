#ifndef DATASTRUCTURES_TARGET_EVENT_H_
#define DATASTRUCTURES_TARGET_EVENT_H_

namespace Scissors {
struct TargetEvent {
  bool special_insertion; // insertions' sequences are known and given by -s
  bool medium_sized_indel;
  bool insertion; // insertions' sequences are unknown

  TargetEvent()
      : special_insertion(false)
      , medium_sized_indel(true)
      , insertion(true)
  {}

  TargetEvent(const bool& s_ins, const bool& indel, const bool& ins)
      : special_insertion(s_ins)
      , medium_sized_indel(indel)
      , insertion(ins)
  {}
};
} // namespace
#endif // DATASTRUCTURES_TARGET_EVENT_H_
