#ifndef DATASTRUCTURES_TARGET_EVENT_H_
#define DATASTRUCTURES_TARGET_EVENT_H_

namespace Scissors {
struct TargetEvent {
  bool special_insertion; // detect the known insertions that are given by -s
  bool special_inversive_insertion; // when special_insertion=true, detect their reverse complement
  bool medium_sized_indel;
  bool insertion; // detect insertions

  TargetEvent()
      : special_insertion(false)
      , special_inversive_insertion(true)
      , medium_sized_indel(true)
      , insertion(true)
  {}

  TargetEvent(const bool& s_ins, const bool& s_inv_ins, const bool& indel, const bool& ins)
      : special_insertion(s_ins)
      , special_inversive_insertion(s_inv_ins)
      , medium_sized_indel(indel)
      , insertion(ins)
  {}
};
} // namespace
#endif // DATASTRUCTURES_TARGET_EVENT_H_
