#ifndef DATASTRUCTURES_THREADDATA_H
#define DATASTRUCTURES_THREADDATA_H

struct ThreadData{
  int id;
  SR_BamInStream* bam_reader;
};

#endif // DATASTRUCTURES_THREADDATA_H
