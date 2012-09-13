#ifndef DATASTRUCTURES_TARGET_REGION_H_
#define DATASTRUCTURES_TARGET_REGION_H_

namespace Scissors {
struct TargetRegion {
  int fragment_length;
  int local_window_size;
  int discovery_window_size;

  TargetRegion()
      : fragment_length(0)
      , local_window_size(1000)
      , discovery_window_size(10000)
  {}
};
} // namespace
#endif // DATASTRUCTURES_TARGET_EVENT_H_
