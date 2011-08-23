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
#include "dataStructures/SR_QueryRegion.h"

//===============================
// Type and constant definition
//===============================

typedef struct HashRegionTable
{
    unsigned int searchBegin;              // lower limit in searching the prevHashArray

    HashRegionArray* pPrevRegions;         // an array hold the hash regions that start at 1bp before

    HashRegionArray* pCurrRegions;         // an array hold the hash regions that start at current bp

    BestRegionArray* pBestCloseRegions;    // an array hold the best hash regions within the closer search region

    BestRegionArray* pBestFarRegions;      // an array hold the best hash regions within the further search region

}HashRegionTable;


//===============================
// Constructors and Destructors
//===============================

HashRegionTable* HashRegionTableAlloc(void);

void HashRegionTableFree(HashRegionTable* pRegionTable);


//===============================
// Interface functions
//===============================

//==================================================================
// function:
//      for each query find the best hash regions in the reference
//
// args:
//      1. pRegionTable: a pointer to a hash region table
//      2. pHashTable: a pointer to a reference hash table
//      3. pQueryRegion: a pointer to a query region
//
// discussion:
//      the best hash region start at each position of the query
//      will be stored at the 'pBestCloseRegions' and the
//      'pBestFarRegions' for close query region and far query
//      region respectively after processing
//==================================================================
void HashRegionTableLoad(HashRegionTable* pRegionTable, const SR_InHashTable* pHashTable, const SR_QueryRegion* pQueryRegion);

//==========================================================
// function:
//      initialize the hash region table for a new query
//
// args:
//      1. pRegionTable: a pointer to a hash region table
//      2. queryLen: length of the query
//==========================================================
void HashRegionTableInit(HashRegionTable* pRegionTable, uint32_t queryLen);

//===========================================================
// function:
//      index the best hash regions with their end position
//      instead of start position
//
// args:
//      1. pRegionTable: a pointer to a hash region table
//
// discussion:
//      by default, the query begin of each best hash region 
//      is its index in the array. for example, in the array 
//      of "pBestCloseRegions", the first element with index 
//      "0" stores the best hash region that starts at query 
//      position "0". This function will reverse each best 
//      hash region so that it stores the best hash region 
//      end instead of begin
//===========================================================
void HashRegionTableReverseBest(HashRegionTable* pRegionTable);

#endif  /*HASHREGIONTABLE_H*/
