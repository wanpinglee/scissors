#ifndef UTILITIES_HASHTABLE_SPECIAL_HASHER_H_
#define UTILITIES_HASHTABLE_SPECIAL_HASHER_H_

#include <string>
extern "C" {
#include "SR_InHashTable.h"
#include "SR_Reference.h"
}

class SpecialHasher {
 public:
  SpecialHasher();
  SpecialHasher(const char* fasta, const int& hash_size = 7);
  ~SpecialHasher();

  // @Description: Setting fasta filename.
  //               Notice that before Load(), the fasta filename should be set.
  //               If the filename is already given in the constructor,
  //               then you don't have to use this function.
  void SetFastaName(const char* fasta) {fasta_ = fasta;};

  // @Description: Setting hash size.
  //               Notice that 1) Default hash size is 7; 
  //               2) before Load(), the hash size should be set.
  //               The size also can be given in the constructor.
  void SetHashSize(const int& hash_size) {hash_size_ = hash_size;};

  // @Description: Loading special references from the fasta file 
  //               and hashing them.
  bool Load();
 private:
  std::string fasta_;
  SR_RefHeader* reference_header_;
  SR_Reference* references_;
  SR_InHashTable* hash_table_;
  int hash_size_;

  void Init();
  SpecialHasher (const SpecialHasher&);
  SpecialHasher& operator= (const SpecialHasher&);
};

#endif //UTILITIES_HASHTABLE_SPECIAL_HASHER_H_
