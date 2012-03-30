#ifndef UTILITIES_HASHTABLE_SPECIAL_HASHER_H_
#define UTILITIES_HASHTABLE_SPECIAL_HASHER_H_

#include <string>
extern "C" {
#include "SR_InHashTable.h"
#include "SR_Reference.h"
}

class SpecialHasher {
 public:
  SpecialHasher(void);
  SpecialHasher(const char* fasta, const int& hash_size = 7);
  ~SpecialHasher(void);

  // @function: Setting fasta filename.
  //            Notice that before Load(), the fasta filename should be set.
  //            If the filename is already given in the constructor,
  //            then you don't have to use this function.
  // @param:    fasta: fasta filename
  void SetFastaName(const char* fasta) {fasta_ = fasta;};

  // @function: Setting hash size.
  //            Notice that 1) Default hash size is 7; 
  //            2) before Load(), the hash size should be set.
  //            The size also can be given in the constructor.
  void SetHashSize(const int& hash_size) {hash_size_ = hash_size;};

  // @function: Loading special references from the fasta file 
  //            and hashing them.
  bool Load(void);

  const SR_Reference* GetReference(void) const {return(is_loaded ? references_ : NULL);};
  const SR_RefHeader* GetReferenceHeader(void) const {return(is_loaded ? reference_header_ : NULL);};
  const SR_InHashTable* GetHashTable(void) const {return(is_loaded ? hash_table_ : NULL);};

 private:
  std::string fasta_;
  SR_RefHeader* reference_header_;
  SR_Reference* references_;
  SR_InHashTable* hash_table_;
  int hash_size_;
  bool is_loaded;

  void Init(void);
  SpecialHasher (const SpecialHasher&);
  SpecialHasher& operator= (const SpecialHasher&);
};

#endif //UTILITIES_HASHTABLE_SPECIAL_HASHER_H_
