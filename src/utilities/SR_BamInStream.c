/*
 * =====================================================================================
 *
 *       Filename:  SR_BamInStream.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/18/2011 05:41:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <assert.h>
#include <limits.h>

#include "samtools/khash.h"
#include "hashTable/common/SR_Error.h"
#include "hashTable/common/SR_Utilities.h"
#include "SR_BamInStream.h"


//===============================
// Type and constant definition
//===============================

// value indicate that we did not load any read
// into the bam array
#define NO_QUERY_YET (-2)

// default capacity of a bam array
#define DEFAULT_BAM_ARRAY_CAP 30

// a mask used to filter out those unwanted reads for split alignments
// it includes proper paired reads, secondar reads, qc-failed reads and duplicated reads
#define SR_BAM_FMASK (BAM_FPROPER_PAIR | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

// initialize a string-to-bam hash table used to retrieve a pair of read
KHASH_MAP_INIT_STR(queryName, bam1_t*);

// array of bam alignments
typedef struct SR_BamArray
{
    bam1_t* data;            // a pointer to the read-in bam alignments

    unsigned int size;       // how many alignments are in the array

    unsigned int capacity;   // maximum number of alignments the array can hold

}SR_BamArray;

// private data structure that holds all bam-input-related 
// information
struct SR_BamInStreamPrvt
{
    bamFile fpBamInput;               // file pointer to a input bam file

    bam_index_t* pBamIndex;           // file pointer to a input bam index file

    bam_iter_t pBamIter;              // pointer to a bam iterator (used for jump)

    SR_BamArray* pBamArray;           // a bam array for all incoming alignments in the bam file

    khash_t(queryName)* nameHash;     // a hash table that holds queryName-bam1_t* pair

    int32_t currRefID;                // the reference ID of the current read-in alignment

    int32_t currPos;                  // the aligned position of the current read-in alignment
};


//===================
// Static functions
//===================

static SR_BamArray* SR_BamArrayAlloc(unsigned int capacity)
{
    SR_BamArray* pBamArray = (SR_BamArray*) malloc(sizeof(SR_BamArray));
    if (pBamArray == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam array object");

    pBamArray->data = (bam1_t*) calloc(capacity, sizeof(bam1_t));
    if (pBamArray->data == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storeage of bam alignments in a bam array object");

    pBamArray->size = 0;
    pBamArray->capacity = capacity;

    return pBamArray;
}

static void SR_BamArrayFree(SR_BamArray* pBamArray)
{
    if (pBamArray != NULL)
    {
        if (pBamArray->data != NULL)
        {
            for (unsigned int i = 0; i != pBamArray->capacity; ++i)
                free(pBamArray->data[i].data);

            free(pBamArray->data);
        }

        free(pBamArray);
    }
}

static void SR_BamArrayStartOver(SR_BamArray* pBamArray)
{
    if (SR_ARRAY_GET_SIZE(pBamArray) <= 1)
        return;

    bam_copy1(SR_ARRAY_GET_PT(pBamArray, 0), SR_ARRAY_GET_LAST_PT(pBamArray));
    pBamArray->size = 1;
}

static int SR_BamInStreamLoadNext(SR_BamInStream* pBamInStream)
{
    if (pBamInStream->pBamArray->size == pBamInStream->pBamArray->capacity)
    {
        pBamInStream->pBamArray->capacity *= 2;
        pBamInStream->pBamArray->data = (bam1_t*) realloc(pBamInStream->pBamArray->data, sizeof(bam1_t) * pBamInStream->pBamArray->capacity);
    }

    int ret = bam_read1(pBamInStream->fpBamInput, SR_ARRAY_GET_PT(pBamInStream->pBamArray, pBamInStream->pBamArray->size));
    ++(pBamInStream->pBamArray->size);

    return ret;
}

static void SR_BamInStreamClear(SR_BamInStream* pBamInStream)
{
    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currPos = NO_QUERY_YET;

    SR_ARRAY_RESET(pBamInStream->pBamArray);
    kh_clear(queryName, pBamInStream->nameHash);
}


//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(const char* bamFilename)
{
    SR_BamInStream* pBamInStream = (SR_BamInStream*) calloc(1, sizeof(struct SR_BamInStreamPrvt));
    if (pBamInStream == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam aux object.");

    pBamInStream->nameHash = kh_init(queryName);
    pBamInStream->pBamArray = SR_BamArrayAlloc(DEFAULT_BAM_ARRAY_CAP);

    pBamInStream->pBamIter = NULL;
    pBamInStream->fpBamInput = bam_open(bamFilename, "r");
    if (pBamInStream->fpBamInput == NULL)
        SR_ErrQuit("ERROR: Cannot open bam file %s for reading.\n", bamFilename);

    pBamInStream->pBamIndex = bam_index_load(bamFilename);
    if (pBamInStream->pBamIndex == NULL)
        SR_ErrMsg("WARNING: Cannot open bam index file for reading. No jump allowed.\n");

    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currPos = NO_QUERY_YET;

    return pBamInStream;
}

void SR_BamInStreamFree(SR_BamInStream* pBamInStream)
{
    if (pBamInStream != NULL)
    {
        kh_destroy(queryName, pBamInStream->nameHash);
        SR_BamArrayFree(pBamInStream->pBamArray);

        bam_iter_destroy(pBamInStream->pBamIter);
        bam_close(pBamInStream->fpBamInput);
        bam_index_destroy(pBamInStream->pBamIndex);

        free(pBamInStream);
    }
}

//======================
// Interface functions
//======================

// jump to a certain chromosome in a bam file
SR_Status SR_BamInStreamJump(SR_BamInStream* pBamInStream, int32_t refID)
{
    // if we do not have the index file return error
    if (pBamInStream->pBamIndex == NULL)
        return SR_ERR;

    // clear the bam array before jump
    SR_BamInStreamClear(pBamInStream);

    // jump and read the first alignment in the given chromosome
    int ret;
    bam_iter_t pBamIter = bam_iter_query(pBamInStream->pBamIndex, refID, 0, INT_MAX);
    ret = bam_iter_read(pBamInStream->fpBamInput, pBamIter, SR_ARRAY_GET_PT(pBamInStream->pBamArray, 0));
    ++(pBamInStream->pBamArray->size);

    bam_iter_destroy(pBamIter);

    // see if we jump to the desired chromosome
    if (ret > 0 && SR_ARRAY_GET(pBamInStream->pBamArray, 0).core.tid == refID)
    {
        pBamInStream->currRefID = refID;
        pBamInStream->currPos = SR_ARRAY_GET(pBamInStream->pBamArray, 0).core.pos;

        return SR_OK;
    }
    else if (ret == -1)
    {
        return SR_OUT_OF_RANGE;
    }
    else
    {
        return SR_ERR;
    }
}

// read the header of a bam file
bam_header_t* SR_BamInStreamReadHeader(SR_BamInStream* pBamInStream)
{
    return bam_header_read(pBamInStream->fpBamInput);
}

// read an alignment from a bam file
SR_Status SR_BamInStreamRead(bam1_t* pAlignment, SR_BamInStream* pBamInStream)
{
    int ret = bam_read1(pBamInStream->fpBamInput, pAlignment);

    if (ret > 0)
        return SR_OK;
    else if (ret == -1)
        return SR_EOF;
    else
        return SR_ERR;
}

// load a unique-orphan pair from a bam file
SR_Status SR_BamInStreamGetPair(bam1_t** ppAnchor, bam1_t** ppOrphan, SR_BamInStream* pBamInStream)
{
    int ret = 1;

    while((ret > 0) && (ret = SR_BamInStreamLoadNext(pBamInStream)) > 0)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        if ((SR_ARRAY_GET_LAST(pBamInStream->pBamArray).core.flag & BAM_FPAIRED) == 0
            || (SR_ARRAY_GET_LAST(pBamInStream->pBamArray).core.flag & SR_BAM_FMASK) != 0
            || strcmp(bam1_qname(SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArray)), "*") == 0)
        {
            SR_ARRAY_POP(pBamInStream->pBamArray);
            continue;
        }

        // update the current ref ID or position if the incoming alignment has a 
        // different value. The name hash and the bam array will be reset
        if (SR_ARRAY_GET_LAST(pBamInStream->pBamArray).core.tid != pBamInStream->currRefID 
            || SR_ARRAY_GET_LAST(pBamInStream->pBamArray).core.pos != pBamInStream->currPos)
        {
            // if we got a different reference ID, we should quit
            if (SR_ARRAY_GET_LAST(pBamInStream->pBamArray).core.tid != pBamInStream->currRefID)
            {
                ret = SR_OUT_OF_RANGE;
            }

            pBamInStream->currRefID = SR_ARRAY_GET_LAST(pBamInStream->pBamArray).core.tid;
            pBamInStream->currPos   = SR_ARRAY_GET_LAST(pBamInStream->pBamArray).core.pos;

            kh_clear(queryName, pBamInStream->nameHash);
            SR_BamArrayStartOver(pBamInStream->pBamArray);

        }

        int khRet = 0;
        khiter_t khIter = kh_put(queryName, pBamInStream->nameHash, bam1_qname(SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArray)), &khRet);

        if (khRet == 0) // we found a pair of alignments 
        {
            // retrieve the unique-orphan pair
            bam1_t* pAlignment = SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArray);
            if ((pAlignment->core.flag & BAM_FUNMAP) == 0)
            {
                (*ppAnchor) = pAlignment;
                (*ppOrphan) = kh_value(pBamInStream->nameHash, khIter);
            }
            else
            {
                (*ppOrphan) = pAlignment;
                (*ppAnchor) = kh_value(pBamInStream->nameHash, khIter);
            }

            // the anchor mate must be aligned uniquely and its mapping quality does not equal to zero
            // the orphan mate must be an unaligned read
            if ((*ppAnchor)->core.qual != 0
                && ((*ppAnchor)->core.flag & BAM_FUNMAP) == 0
                && ((*ppOrphan)->core.flag & BAM_FUNMAP) != 0)
            {
                ret = SR_OK;
            }
        }
        else // not finding corresponding mate, save the current value and move on
        {
            kh_value(pBamInStream->nameHash, khIter) = SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArray);
        }
    }

    if (ret < 0)
    {
        if (ret != SR_EOF && ret != SR_OUT_OF_RANGE)
            return SR_ERR;
    }

    return ret;
}
