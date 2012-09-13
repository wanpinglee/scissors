#ifndef DATASTRUCTURES_ZA_TAG_H_
#define DATASTRUCTURES_ZA_TAG_H_

enum ZaSvType {
  kZaSvTypeNone    = 0, 
  kZaSvTypeDel     = 1, 
  kZaSvTypeInv     = 2, 
  kZaSvTypeSpecial = 3
};

static char * ZaSvString[] = {
  "NON",
  "DEL",
  "INV",
  "SPE"
};

struct Zatag {
  bool at; // true: '@'; false: '&'
  int mq;
  ZaSvType type;
  char special_letter[2];
  int mappings;
  string cigar;
  string md;

  Zatag()
      : at(false)
      , mq(0)
      , type(kZaSvTypeNone)
      , special_letter("\0")
      , mappings(0)
      , cigar()
      , md()
  {};

  Clear() {
    at = false;
    mq = 0;
    type = kZaSvTypeNone;
    special_letter = "\0";
    mappings = 0;
    cigar.clear();
    md = clear();
  }
}

#endif // DATASTRUCTURES_ZA_TAG_H_
