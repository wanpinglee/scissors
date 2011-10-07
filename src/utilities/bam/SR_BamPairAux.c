/*
 * =====================================================================================
 *
 *       Filename:  SR_BamPairAux.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/17/2011 09:43:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "outsources/samtools/bam.h"
#include "SR_BamPairAux.h"


// default soft clipping tolerance
#define DEFAULT_SC_TOLERANCE 0.2

// alignment status
enum AlignmentStatus
{
    NEITHER_GOOD = -1,    // neither a good anchor nor a good orphan candidate

    GOOD_ANCHOR  = 0,     // a good anchor candidate

    GOOD_ORPHAN  = 1      // a good orphan candidate
};

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

SR_Bool SR_IsUniqueOrphanPair(SR_BamNode** ppAnchor, SR_BamNode** ppOrphan, double scTolerance)
{
    int anchorStatus = SR_CheckSC(&((*ppAnchor)->alignment), scTolerance);
    int orphanStatus = SR_CheckSC(&((*ppOrphan)->alignment), scTolerance);

    if ((anchorStatus == NEITHER_GOOD || orphanStatus == NEITHER_GOOD)
        || (anchorStatus == GOOD_ANCHOR && orphanStatus == GOOD_ANCHOR)
        || (anchorStatus == GOOD_ORPHAN && orphanStatus == GOOD_ORPHAN))
    {
        return FALSE;
    }
    else if (anchorStatus == GOOD_ORPHAN)
        SR_SWAP(*ppAnchor, *ppOrphan, SR_BamNode*);

    return TRUE;
}

SR_Status SR_LoadUniquOrphanPairs(SR_BamInStream* pBamInStream, unsigned int threadID, double scTolerance)
{
    SR_BamNode* pAnchor = NULL;
    SR_BamNode* pOrphan = NULL;

    SR_Status readerStatus = SR_OK;
    SR_Status bufferStatus = SR_OK;
    while ((readerStatus = SR_BamInStreamLoadPair(&pAnchor, &pOrphan, NULL, pBamInStream)) == SR_OK)
    {
        if (SR_IsUniqueOrphanPair(&pAnchor, &pOrphan, scTolerance))
        {
            SR_BamInStreamPush(pBamInStream, pAnchor, threadID);
            bufferStatus = SR_BamInStreamPush(pBamInStream, pOrphan, threadID);
        }
        else
        {
            SR_BamInStreamRecycle(pBamInStream, pAnchor);
            SR_BamInStreamRecycle(pBamInStream, pOrphan);
        }

        if (bufferStatus == SR_FULL)
            break;
    }

    return readerStatus;
}

SR_Status SR_BamPairStatsLoad(SR_BamPairStats* pPairStats, SR_BamNode* pUpAlgn, SR_BamNode* pDownAlgn)
{
    SR_SingleMode upMode = SR_BamGetSingleMode(&(pUpAlgn->alignment));
    SR_SingleMode downMode = SR_BamGetSingleMode(&(pDownAlgn->alignment));
    pPairStats->pairMode = SR_BamGetPairMode(upMode, downMode);
    if (pPairStats->pairMode < 0)
        return SR_ERR;

    pPairStats->fragLen = abs(pUpAlgn->alignment.core.isize);

    static const char tagRG[2] = {'R', 'G'};
    uint8_t* rgPos = bam_aux_get(&(pUpAlgn->alignment), tagRG);
    if (rgPos != NULL)
        pPairStats->RG = bam_aux2Z(rgPos);
    else
        pPairStats->RG = NULL;

    return SR_OK;
}
