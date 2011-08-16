
#ifndef SRC_UTILITIES_BAMUTILITIES_H_
#define SRC_UTILITIES_BAMUTILITIES_H_

#include <string>

using namespace std;

// Given sequence, generate bam-format encoded sequence 
// and store it in encodedSequence
void EncodeQuerySequence( string& encodedSequence, const string& sequence );


#endif // SRC_UTILITIES_BAMUTILITIES_H_
