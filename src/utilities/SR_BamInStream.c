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

#include <limits.h>
#include <assert.h>

#include "samtools/khash.h"
#include "hashTable/common/SR_Error.h"
#include "SR_BamInStream.h"
#include "hashTable/common/SR_Utilities.h"


//===============================
// Type and constant definition
//===============================

// value indicate that we did not load any read
// into the bam array
#define NO_QUERY_YET (-2)

// default capacity of a bam array
#define DEFAULT_BAM_ARRAY_CAP 50

// default soft clipping tolerance
#define DEFAULT_SC_TOLERANCE 0.2

// alignment status
enum AlignmentStatus
{
    NEITHER_GOOD = -1,    // neither a good anchor nor a good orphan candidate

    GOOD_ANCHOR  = 0,     // a good anchor candidate

    GOOD_ORPHAN  = 1      // a good orphan candidate
};

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

// private data structure that holds all bam-input-related information
struct SR_BamInStreamPrvt
{
    bamFile fpBamInput;                   // file pointer to a input bam file

    bam_index_t* pBamIndex;               // file pointer to a input bam index file

    SR_BamArray* pBamArrayPrev;           // a bam array for all incoming alignments in the previous bin

    SR_BamArray* pBamArrayCurr;           // a bam array for all incoming alignments in the current bin

    khash_t(queryName)* nameHashPrev;     // a hash table that holds queryName-bam1_t* pair in previous bin

    khash_t(queryName)* nameHashCurr;     // a hash table that holds queryName-bam1_t* pair in current bin

    int32_t currRefID;                    // the reference ID of the current read-in alignment

    int32_t currBinPos;                   // the start position of current bin (0-based)

    uint32_t binLen;                      // the length of bin

    double scTolerance;                   // soft clipping tolerance rate. 
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
    if (pBamInStream->pBamArrayCurr->size == pBamInStream->pBamArrayCurr->capacity)
    {
        pBamInStream->pBamArrayCurr->capacity *= 2;
        pBamInStream->pBamArrayCurr->data = (bam1_t*) realloc(pBamInStream->pBamArrayCurr->data, sizeof(bam1_t) * pBamInStream->pBamArrayCurr->capacity);
    }

    ++(pBamInStream->pBamArrayCurr->size);
    int ret = bam_read1(pBamInStream->fpBamInput, SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr));

    return ret;
}

static int SR_CheckSC(bam1_t* pAlignment, double scTolerance)
{
    if ((pAlignment->core.flag & BAM_FUNMAP) != 0)
        return GOOD_ORPHAN;

    unsigned int scLimit = scTolerance * pAlignment->core.l_qseq;
    uint32_t* cigar = bam1_cigar(pAlignment);

    SR_Bool isHeadSC = FALSE;
    SR_Bool isTailSC = FALSE;

    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP 
        && (cigar[0] >> BAM_CIGAR_SHIFT) >= scLimit)
    {
        isHeadSC = TRUE;
    }

    unsigned int lastIndex = pAlignment->core.n_cigar - 1;
    if ((cigar[lastIndex] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP 
        && (cigar[lastIndex] >> BAM_CIGAR_SHIFT) >= scLimit)
    {
        isTailSC = TRUE;
    }

    if (!isHeadSC && !isTailSC) 
    {
        if (pAlignment->core.qual != 0)
            return GOOD_ANCHOR;
        else
            return NEITHER_GOOD;
    }
    else if (isHeadSC && isTailSC)
        return NEITHER_GOOD;
    else
        return GOOD_ORPHAN;

}

static SR_Bool SR_IsQualifiedPair(bam1_t** ppAnchor, bam1_t** ppOrphan, double scTolerance)
{
    int anchorStatus = SR_CheckSC((*ppAnchor), scTolerance);
    int orphanStatus = SR_CheckSC((*ppOrphan), scTolerance);

    if ((anchorStatus == NEITHER_GOOD || orphanStatus == NEITHER_GOOD)
        || (anchorStatus == GOOD_ANCHOR && orphanStatus == GOOD_ANCHOR)
        || (anchorStatus == GOOD_ORPHAN && orphanStatus == GOOD_ORPHAN))
    {
        return FALSE;
    }
    else if (anchorStatus == GOOD_ORPHAN)
        SR_SWAP(*ppAnchor, *ppOrphan, bam1_t*);

    return TRUE;
}

static void SR_BamInStreamClear(SR_BamInStream* pBamInStream)
{
    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->currRefID = NO_QUERY_YET;

    kh_clear(queryName, pBamInStream->nameHashPrev);
    kh_clear(queryName, pBamInStream->nameHashCurr);

    SR_ARRAY_RESET(pBamInStream->pBamArrayPrev);
    SR_ARRAY_RESET(pBamInStream->pBamArrayCurr);
}


//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(const char* bamFilename, uint32_t binLen, double scTolerance)
{
    SR_BamInStream* pBamInStream = (SR_BamInStream*) calloc(1, sizeof(struct SR_BamInStreamPrvt));
    if (pBamInStream == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam aux object.");

    pBamInStream->nameHashPrev = kh_init(queryName);
    pBamInStream->nameHashCurr = kh_init(queryName);

    pBamInStream->pBamArrayPrev = SR_BamArrayAlloc(DEFAULT_BAM_ARRAY_CAP);
    pBamInStream->pBamArrayCurr = SR_BamArrayAlloc(DEFAULT_BAM_ARRAY_CAP);

    pBamInStream->fpBamInput = bam_open(bamFilename, "r");
    if (pBamInStream->fpBamInput == NULL)
        SR_ErrQuit("ERROR: Cannot open bam file %s for reading.\n", bamFilename);

    pBamInStream->pBamIndex = bam_index_load(bamFilename);
    if (pBamInStream->pBamIndex == NULL)
        SR_ErrMsg("WARNING: Cannot open bam index file for reading. No jump allowed.\n");

    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->binLen = binLen;

    pBamInStream->scTolerance = scTolerance;

    return pBamInStream;
}

void SR_BamInStreamFree(SR_BamInStream* pBamInStream)
{
    if (pBamInStream != NULL)
    {
        kh_destroy(queryName, pBamInStream->nameHashPrev);
        kh_destroy(queryName, pBamInStream->nameHashCurr);

        SR_BamArrayFree(pBamInStream->pBamArrayPrev);
        SR_BamArrayFree(pBamInStream->pBamArrayCurr);

        bam_close(pBamInStream->fpBamInput);
        bam_index_destroy(pBamInStream->pBamIndex);

        free(pBamInStream);
    }
}

SR_BamHeader* SR_BamHeaderAlloc(void)
{
    SR_BamHeader* pNewHeader = (SR_BamHeader*) calloc(1, sizeof(SR_BamHeader));
    if (pNewHeader == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam header object");

    return pNewHeader;
}

void SR_BamHeaderFree(SR_BamHeader* pBamHeader)
{
    if (pBamHeader != NULL)
    {
        free(pBamHeader->pMD5s);
        bam_header_destroy(pBamHeader->pOrigHeader);

        free(pBamHeader);
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

    ++(pBamInStream->pBamArrayCurr->size);
    ret = bam_iter_read(pBamInStream->fpBamInput, pBamIter, SR_ARRAY_GET_FIRST_PT(pBamInStream->pBamArrayCurr));

    bam_iter_destroy(pBamIter);

    // see if we jump to the desired chromosome
    if (ret > 0 && SR_ARRAY_GET_FIRST(pBamInStream->pBamArrayCurr).core.tid == refID)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        if ((SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.flag & SR_BAM_FMASK) != 0
            || strcmp(bam1_qname(SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr)), "*") == 0)
        {
            SR_ARRAY_POP(pBamInStream->pBamArrayCurr);
            pBamInStream->currBinPos = NO_QUERY_YET;
        }
        else
        {
            int khRet = 0;
            khiter_t khIter = kh_put(queryName, pBamInStream->nameHashCurr, bam1_qname(SR_ARRAY_GET_FIRST_PT(pBamInStream->pBamArrayCurr)), &khRet);

            if (khRet != 0)
                kh_value(pBamInStream->nameHashCurr, khIter) = SR_ARRAY_GET_FIRST_PT(pBamInStream->pBamArrayCurr);
            else
                return SR_ERR;

            pBamInStream->currBinPos = SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.pos;
        }

        pBamInStream->currRefID = refID;
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
SR_BamHeader* SR_BamInStreamLoadHeader(SR_BamInStream* pBamInStream)
{
    bam_header_t* pOrigHeader = bam_header_read(pBamInStream->fpBamInput);
    if (pOrigHeader == NULL)
        return NULL;

    SR_BamHeader* pBamHeader = SR_BamHeaderAlloc();

    pBamHeader->pOrigHeader = pOrigHeader;

    pBamHeader->pMD5s = (const char**) calloc(pOrigHeader->n_targets, sizeof(char*));
    if (pBamHeader->pMD5s == NULL)
        SR_ErrQuit("ERROR: Not enough memory for md5 string");

    unsigned int numMD5 = 0;
    for (const char* md5Pos = pOrigHeader->text; numMD5 <= pOrigHeader->n_targets && (md5Pos = strstr(md5Pos, "M5:")) != NULL; ++numMD5, ++md5Pos)
    {
        pBamHeader->pMD5s[numMD5] = md5Pos + 3;
    }

    if (numMD5 != pOrigHeader->n_targets)
    {
        free(pBamHeader->pMD5s);
        pBamHeader->pMD5s = NULL;

        if (numMD5 != 0)
            SR_ErrMsg("WARNING: Number of MD5 string is not consistent with number of chromosomes.");
    }

    return pBamHeader;
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

// get the current reference ID
int32_t SR_BamInStreamGetRefID(const SR_BamInStream* pBamInStream)
{
    return (pBamInStream->currRefID);
}

// load a unique-orphan pair from a bam file
SR_Status SR_BamInStreamGetPair(bam1_t** ppAnchor, bam1_t** ppOrphan, SR_BamInStream* pBamInStream)
{
    int ret = 1;

    while((ret > 0) && (ret = SR_BamInStreamLoadNext(pBamInStream)) > 0)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        if ((SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.flag & SR_BAM_FMASK) != 0
            || strcmp(bam1_qname(SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr)), "*") == 0)
        {
            SR_ARRAY_POP(pBamInStream->pBamArrayCurr);
            continue;
        }

        // update the current ref ID or position if the incoming alignment has a 
        // different value. The name hash and the bam array will be reset
        if (SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.tid != pBamInStream->currRefID
            || SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.pos >= pBamInStream->currBinPos + 2 * pBamInStream->binLen)
        {
            pBamInStream->currRefID  = SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.tid;
            pBamInStream->currBinPos = SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.pos;

            kh_clear(queryName, pBamInStream->nameHashPrev);
            kh_clear(queryName, pBamInStream->nameHashCurr);

            SR_ARRAY_RESET(pBamInStream->pBamArrayPrev);
            SR_BamArrayStartOver(pBamInStream->pBamArrayCurr);
        }
        else if (SR_ARRAY_GET_LAST(pBamInStream->pBamArrayCurr).core.pos >= pBamInStream->currBinPos + pBamInStream->binLen)
        {
            pBamInStream->currBinPos += pBamInStream->binLen;

            kh_clear(queryName, pBamInStream->nameHashPrev);
            SR_SWAP(pBamInStream->nameHashPrev, pBamInStream->nameHashCurr, khash_t(queryName)*);

            SR_ARRAY_RESET(pBamInStream->pBamArrayPrev);
            SR_SWAP(pBamInStream->pBamArrayPrev, pBamInStream->pBamArrayCurr, SR_BamArray*);

            bam_copy1(SR_ARRAY_GET_FIRST_PT(pBamInStream->pBamArrayCurr), SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayPrev));
            ++(pBamInStream->pBamArrayCurr->size);

            SR_ARRAY_POP(pBamInStream->pBamArrayPrev);
        }

        khiter_t khIter;
        (*ppAnchor) = NULL;
        (*ppOrphan) = NULL;

        khIter = kh_get(queryName, pBamInStream->nameHashPrev, bam1_qname(SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr)));
        if (khIter != kh_end(pBamInStream->nameHashPrev))
        {
            (*ppAnchor) = kh_value(pBamInStream->nameHashPrev, khIter);
            (*ppOrphan) = SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr);

            if (SR_IsQualifiedPair(ppAnchor, ppOrphan, pBamInStream->scTolerance))
                ret = SR_OK;
        }
        else
        {
            int khRet = 0;
            khIter = kh_put(queryName, pBamInStream->nameHashCurr, bam1_qname(SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr)), &khRet);

            if (khRet == 0) // we found a pair of alignments 
            {
                (*ppAnchor) = kh_value(pBamInStream->nameHashCurr, khIter);
                (*ppOrphan) = SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr);

                if (SR_IsQualifiedPair(ppAnchor, ppOrphan, pBamInStream->scTolerance))
                    ret = SR_OK;
            }
            else // not finding corresponding mate, save the current value and move on
            {
                kh_value(pBamInStream->nameHashCurr, khIter) = SR_ARRAY_GET_LAST_PT(pBamInStream->pBamArrayCurr);
            }
        }
    }

    if (ret < 0)
    {
        if (ret != SR_EOF)
            return SR_ERR;
    }

    return ret;
}
