/*
 * =====================================================================================
 *
 *       Filename:  SR_BamPairAux.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/17/2011 09:43:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_BAMPAIRAUX_H
#define  SR_BAMPAIRAUX_H


#include "utilities/common/SR_Types.h"
#include "SR_BamMemPool.h"
#include "SR_BamInStream.h"
#include "SR_FragLenDstrb.h"
#include "utilities/common/SR_Utilities.h"


//===============================
// Type and constant definition
//===============================

static const unsigned int SR_UNIQUE_ORPHAN_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);

static const unsigned int SR_NORMAL_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FMUNMAP);

//======================
// Interface functions
//======================

static inline SR_Bool SR_CommonFilter(SR_BamNode* pBamNode, const void* filterData)
{
    if ((pBamNode->alignment.core.flag & BAM_FPAIRED) == 0
        || strcmp(bam1_qname(&(pBamNode->alignment)), "*") == 0
        || (pBamNode->alignment.core.flag & SR_UNIQUE_ORPHAN_FMASK) != 0
        || (pBamNode->alignment.core.flag | (BAM_FUNMAP | BAM_FMUNMAP)) == pBamNode->alignment.core.flag)
    {
        return TRUE;
    }
    
    return FALSE;
}

static inline SR_Bool SR_NormalFilter(SR_BamNode* pBamNode, const void* filterData)
{
    if ((pBamNode->alignment.core.flag & BAM_FPAIRED) == 0
        || strcmp(bam1_qname(&(pBamNode->alignment)), "*") == 0
        || (pBamNode->alignment.core.flag & SR_NORMAL_FMASK) != 0
        || (pBamNode->alignment.core.isize == 0))
    {
        return TRUE;
    }

    // any reads aligned to different chromosome will be kept as SV candidates
    if (pBamNode->alignment.core.tid != pBamNode->alignment.core.mtid)
        return TRUE;

    return FALSE;
}

//==================================================================
// function:
//      check if a pair of alignment is qualified unique-orphan 
//      pair
//
// args:
//      ppAnchor: a pointer to a pointer of bam node with anchor
//                alignment
//      ppOrphan: a pointer to a pointer of bam node with orphan
//                alignment
//      scTolerance: soft clipping tolerance
//
// return:
//      if they are qualified unique orphan pair return TRUE
//      else return FALSE
//==================================================================
SR_AlgnType SR_GetAlignmentType(SR_BamNode** ppAlgnOne, SR_BamNode** ppAlgnTwo, double scTolerance, double maxMismatchRate, unsigned char minMQ);

//==================================================================
// function:
//      load a certain number of unique orphan pairs into a buffer
//      associated with a thread
//
// args:
//      1. pBamInStream: a pointer to a bam in stream object
//      2. threadID: the ID of the thread
//      3. scTolerance: soft clipping tolerance
//
// return:
//      status of the bam in stream. if we reach the end of a
//      chromosome, return SR_OUT_OF_RANGE; if we reach the end of
//      the file, return SR_EOF; if an error happens, return
//      SR_ERR; else return SR_OK
//==================================================================
SR_Status SR_LoadAlgnPairs(SR_BamInStream* pBamInStream, SR_FragLenDstrb* pDstrb, unsigned int threadID, double scTolerance, double maxMismatchRate, unsigned char minMQ);

//====================================================================
// function:
//      check if a pair of read is normal
//
// args:
//      1. ppUpAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with smaller coordinate
//      1. ppDownAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with greater coordinate
//
// return:
//      if the pair of read is normal, return TRUE; else, return
//      FALSE
//=====================================================================
static inline SR_Bool SR_IsNormalPair(SR_BamNode** ppUpAlgn, SR_BamNode** ppDownAlgn, unsigned short minMQ)
{
    if (((*ppUpAlgn)->alignment.core.flag & BAM_FUNMAP) != 0
        || ((*ppDownAlgn)->alignment.core.flag & BAM_FUNMAP) != 0
        || (*ppUpAlgn)->alignment.core.qual < minMQ 
        || (*ppDownAlgn)->alignment.core.qual < minMQ)
    {
        return FALSE;
    }

    return TRUE;
}
    
#endif  /*SR_BAMPAIRAUX_H*/
