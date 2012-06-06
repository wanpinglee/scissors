#ifndef DATASTRUCTURES_TARGET_EVENT_H_
#define DATASTRUCTURES_TARGET_EVENT_H_

namespace Scissors {
struct TargetEvent {
  bool special_insertion;
  bool medium_sized_indel;

  TargetEvent()
      : special_insertion(false)
      , medium_sized_indel(true)
  {}

  TargetEvent(const bool& s_insertion, const bool& m_indel)
      : special_insertion(s_insertion)
      , medium_sized_indel(m_indel)
  {}
};
} // namespace
#endif // DATASTRUCTURES_TARGET_EVENT_H_
