/*
 * =====================================================================================
 *
 *       Filename:  HashRegionTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/03/2011 15:22:49
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  HASHREGIONTABLE_H
#define  HASHREGIONTABLE_H

#include "SR_HashRegionArray.h"
#include "SR_InHashTable.h"

//===============================
// Type and constant definition
//===============================

typedef struct HashRegionTable
{
    unsigned int searchBegin;              // lower limit in searching the prevHashArray

    double edgeTolPercent;                 // the best hash regions should start at a read no later than this position

    HashRegionArray* pPrevRegions;         // an array hold the hash regions that start at 1bp before

    HashRegionArray* pCurrRegions;        // an array hold the hash regions that start at current bp

    BestRegionArray* pBestCloseRegions;    // an array hold the best hash regions within the closer search region

    BestRegionArray* pBestFarRegions;      // an array hold the best hash regions within the further search region

}HashRegionTable;

typedef struct QueryRegion
{
    const char* query;              // sequencing read

    unsigned short queryLen;        // length of the read

    uint32_t refBegin;              // the beginning position of search region for the first partial alignment (closer to the anchor mate)

    uint32_t closeRefBound;         // the end position of the search region for the first partial alignment (closer to the anchor mate)

    uint32_t farRefBound;           // the end position of search region for the second partial alignment (further from the anchor mate)

}QueryRegion;


//===============================
// Constructors and Destructors
//===============================

HashRegionTable* HashRegionTableAlloc(double edgeTolPercent);

void HashRegionTableFree(HashRegionTable* pRegionTable);


//===============================
// Non-constant methods
//===============================

// for each query find the best hash regions in the reference
void HashRegionTableLoad(HashRegionTable* pRegionTable, const SR_InHashTable* pHashTable, const QueryRegion* pQueryRegion);

// by default, the query begin of each best hash region is its index in the array
// for example, in the array of "pBestCloseRegions", the first element with index "0" stores the best hash region
// that starts at query position "0". This function will reverse each best hash region so that it stores the best
// hash region end in that position instead of begin
void HashRegionTableReverseBest(HashRegionTable* pRegionTable);

// initialize the hash region table for a new query
void HashRegionTableInit(HashRegionTable* pRegionTable, unsigned short queryLen);

#endif  /*HASHREGIONTABLE_H*/
