
/*
 * =====================================================================================
 *
 *       Filename:  HashRegionArray.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/04/2011 11:00:03
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  HASHARRAY_H
#define  HASHARRAY_H

#include <stdint.h>

#include "hashTable/common/SR_Types.h"

//===============================
// Type and constant definition
//===============================

#define MAX_BEST_REF_BEGINS 3

// default capacity of a hash region array
static const uint32_t DEFAULT_HASH_ARR_CAPACITY = 50;

// structure hold the information of a hash region
typedef struct HashRegion
{
    uint32_t refBegin;     // begin position of a hash region at the reference

    uint32_t  queryBegin;   // begin position of a hash region at a read

    uint32_t  length;       // length of a hash region

}HashRegion;

// structure to hold the information about the best hash region start at a certain position in a read
typedef struct BestRegion
{
    uint32_t refBegins[MAX_BEST_REF_BEGINS];

    uint32_t queryBegin;

    uint32_t length;

    uint32_t numPos;

}BestRegion;

// hash region array
typedef struct HashRegionArray
{
    HashRegion* data;      // array of hash regions 

    unsigned int size;     // size of the array

    unsigned int capacity; // maximum number of objects that can be held in the array

}HashRegionArray;

typedef struct BestRegionArray
{
    BestRegion* data;

    unsigned int size;

    unsigned int capacity;

}BestRegionArray;


#endif  /*HASHARRAY_H*/
