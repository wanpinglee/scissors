/*
 * =====================================================================================
 *
 *       Filename:  SR_InHashTable.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/13/2011 00:49:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>

#include "hashTable/common/SR_Error.h"
#include "SR_InHashTable.h"

SR_InHashTable* SR_InHashTableAlloc(unsigned char hashSize)
{
    SR_InHashTable* pNewTable = (SR_InHashTable*) malloc(sizeof(SR_InHashTable));
    if (pNewTable == NULL)
        SR_ErrSys("ERROR: Not enough memory for a reference hash table object.\n");
    
    pNewTable->chr = 0;
    pNewTable->hashSize = hashSize;
    pNewTable->md5[MD5_STR_LEN] = '\0';

    pNewTable->highEndMask = GET_HIGH_END_MASK(hashSize);
    pNewTable->numHashes = (uint32_t) 1 << (2 * hashSize);
    pNewTable->indices = (uint32_t*) malloc(sizeof(uint32_t) * pNewTable->numHashes);
    if (pNewTable->indices == NULL)
        SR_ErrSys("ERROR: Not enough memory for the hash index array in a hash table object.\n");

    pNewTable->numPos = 0;
    pNewTable->hashPos = NULL;

    return pNewTable;
}

void SR_InHashTableFree(SR_InHashTable* pHashTable)
{
    if (pHashTable != NULL)
    {
        free(pHashTable->hashPos);
        free(pHashTable->indices);

        free(pHashTable);
    }
}

SR_Bool SR_InHashTableRead(SR_InHashTable* pHashTable, FILE* hashTableInput)
{
    size_t readSize = 0;

    readSize = fread(&(pHashTable->chr), sizeof(unsigned char), 1, hashTableInput);
    if (readSize != 1)
    {
        if (feof(hashTableInput))
            return FALSE;
        else
            SR_ErrSys("ERROR: Cannot read the chromsome number from the hash table file.\n");
    }

    readSize = fread(pHashTable->md5, sizeof(char), MD5_STR_LEN, hashTableInput);
    if (readSize != MD5_STR_LEN)
        SR_ErrSys("ERROR: Cannot read the md5 checksum string from the hash table file.\n");

    readSize = fread(pHashTable->indices, sizeof(uint32_t), pHashTable->numHashes, hashTableInput);
    if (readSize != pHashTable->numHashes)
        SR_ErrSys("ERROR: Cannot read the indices from the hash table file.\n");

    readSize = fread(&(pHashTable->numPos), sizeof(uint32_t), 1, hashTableInput);
    if (readSize != 1)
        SR_ErrSys("ERROR: Cannot read the total number of hash positions from the hash table file.\n");

    free(pHashTable->hashPos);
    pHashTable->hashPos = (uint32_t*) malloc(sizeof(uint32_t) * pHashTable->numPos);

    readSize = fread(pHashTable->hashPos, sizeof(uint32_t), pHashTable->numPos, hashTableInput);
    if (readSize != pHashTable->numPos)
        SR_ErrSys("ERROR: Cannot read the hash positions from the hash table.\n");

    return TRUE;
}


SR_Bool SR_InHashTableSearch(HashPosView* hashPosView, const SR_InHashTable* pHashTable, uint32_t hashKey)
{
    if(hashKey >= pHashTable->numHashes)
        SR_ErrSys("ERROR: Invalid hash key.\n");

    uint32_t index = pHashTable->indices[hashKey];
    uint32_t nextIndex = hashKey == (pHashTable->numHashes - 1) ? pHashTable->numPos : pHashTable->indices[hashKey + 1];

    if (index == nextIndex)
        return FALSE;
    
    hashPosView->size = nextIndex - index;
    hashPosView->data = pHashTable->hashPos + index;

    return TRUE;
}
