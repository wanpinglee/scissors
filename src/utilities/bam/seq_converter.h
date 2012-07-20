#ifndef UTILITIES_BAM_SEQ_CONVERTER_H_
#define UTILITIES_BAM_SEQ_CONVERTER_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
void GetReverseComplementSequence(const uint8_t* original, 
                                  const int length, 
              	 		  uint8_t* reverse);

void GetComplementSequence(const uint8_t* original,
                           const int length,
			   uint8_t* reverse);

void GetInverseSequence(const uint8_t* original,
                        const int length,
			uint8_t* reverse);
#ifdef __cplusplus
}
#endif
#endif
