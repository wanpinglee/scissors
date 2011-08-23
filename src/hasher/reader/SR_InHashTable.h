/*
 * =====================================================================================
 *
 *       Filename:  SR_InHashTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/12/2011 22:32:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  SR_INHASHTABLE_H
#define  SR_INHASHTABLE_H

#include <stdio.h>
#include <stdint.h>

#include "hasher/common/SR_Types.h"


//===============================
// Type and constant definition
//===============================

// generate a mask to clear the highest 2 bits in a hash key (the leftmost base pair)
#define GET_HIGH_END_MASK(hashSize) ((uint32_t) 0xffffffff >> (34 - (2 * (hashSize))))

#define READ_HASH_SIZE(hashSize, hashTableInput)                                           \
    do                                                                                     \
    {                                                                                      \
        size_t readSize = 0;                                                               \
        readSize = fread(&(hashSize), sizeof(unsigned char), 1, hashTableInput);           \
        if (readSize != 1)                                                                 \
            SR_ErrQuit("ERROR: Cannot read hash size from the hash table file.\n ");       \
    }while(0)

typedef struct HashPosView
{
    const uint32_t* data;     // the address where a certain hash start in the "hashPos" array in the "SR_InHashTable" object 
    unsigned int size;        // total number of hash position found in the reference for a certain hash 

}HashPosView;

// input format of reference hash table
typedef struct SR_InHashTable
{
    unsigned char chr;             // chromosome of the current reference hash table

    unsigned char hashSize;        // size of hash

    char md5[MD5_STR_LEN + 1];     // md5 checksum string

    uint32_t* hashPos;             // positions of hashes found in the reference sequence

    uint32_t* indices;             // index of a given hash in the "hashPos" array

    uint32_t  highEndMask;         // a mask to clar the highest 2 bits in a hash key

    uint32_t  numPos;              // total number of hash positions found in reference

    uint32_t  numHashes;           // total number of different hashes

}SR_InHashTable;


//===============================
// Constructors and Destructors
//===============================

SR_InHashTable* SR_InHashTableAlloc(unsigned char hashSize);

void SR_InHashTableFree(SR_InHashTable* pHashTable);


//===============================
// Non-constant methods
//===============================

SR_Bool SR_InHashTableRead(SR_InHashTable* pHashTable, FILE* hashTableInput);

SR_Bool SR_InHashTableSearch(HashPosView* hashPosView, const SR_InHashTable* pHashTable, uint32_t hashKey);


#endif  /*SR_INHASHTABLE_H*/
