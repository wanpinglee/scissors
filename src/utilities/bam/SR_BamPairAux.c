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

#define DEFAULT_MAX_MISMATCH_RATE 0.1

#define BAM_CMISMATCH 8

// alignment status
enum AlignmentStatus
{
    NONE_GOOD      = -1,    // neither a good anchor nor a good orphan candidate
    
    GOOD_ANCHOR    = 0,     // a good anchor candidate
    
    GOOD_ORPHAN    = 1,     // a good orphan candidate
    
    GOOD_SOFT      = 2,     // a good soft clipping candidate
    
    GOOD_MULTIPLE  = 3,     // a good multiple aligned candidate

    GOOD_POOR      = 4,     // a good candidate that poorly mapped
};

static double SR_GetMismatchRate(const bam1_t* pAlignment)
{
    uint32_t* cigar = bam1_cigar(pAlignment);
    unsigned int numMismatch = 0;
    for (unsigned i = 0; i != pAlignment->core.n_cigar; ++i)
    {
        int type = (cigar[i] & BAM_CIGAR_MASK);
        //if (type == BAM_CINS || type == BAM_CDEL || type == BAM_CSOFT_CLIP || type == BAM_CMISMATCH)
	if (type == BAM_CINS || type == BAM_CDEL || type == BAM_CMISMATCH)
        {
            numMismatch += (cigar[i] >> BAM_CIGAR_SHIFT);
        }
    }

    return ((double) numMismatch / pAlignment->core.l_qseq);
}

static int SR_CheckAlignment(const bam1_t* pAlignment, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    if ((pAlignment->core.flag & BAM_FUNMAP) != 0) // the read itself is unmapped
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

    // the alignment passes the soft-clip checker
    if (!isHeadSC && !isTailSC) 
    {
        double mismatchRate = SR_GetMismatchRate(pAlignment);

        #ifdef VERBOSE_DEBUG
	  if (pAlignment->core.qual < minMQ) 
	    fprintf(stderr, "\tQuality low: %d < %d\n", pAlignment->core.qual, minMQ);
	  if (mismatchRate > maxMismatchRate) 
	    fprintf(stderr, "\tToo many X: %1.2f > %1.2f\n", mismatchRate, maxMismatchRate);
	#endif

	if ((pAlignment->core.qual >= minMQ) && (mismatchRate <= maxMismatchRate))
            return GOOD_ANCHOR;
	else
	    return GOOD_POOR; // Low quality alignment
    }
    // the alignment has too many soft clips
    else {
        return GOOD_SOFT;
    }
}

SR_AlgnType SR_GetAlignmentType(SR_BamNode** ppAlgnOne, SR_BamNode** ppAlgnTwo, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    #ifdef VERBOSE_DEBUG
      fprintf(stderr,"%s\n", bam1_qname(&((*ppAlgnOne)->alignment)));
    #endif

    int firstType = SR_CheckAlignment(&((*ppAlgnOne)->alignment), scTolerance, maxMismatchRate, minMQ);
    
    #ifdef VERBOSE_DEBUG
      fprintf(stderr,"mate1 type: %d\n", firstType);
    #endif
    
    int secondType = SR_CheckAlignment(&((*ppAlgnTwo)->alignment), scTolerance, maxMismatchRate, minMQ);
    
    #ifdef VERBOSE_DEBUG
      fprintf(stderr,"mate2 type: %d\n", secondType);
    #endif

    if (firstType != GOOD_ANCHOR && secondType == GOOD_ANCHOR)
    {
        SR_SWAP(*ppAlgnOne, *ppAlgnTwo, SR_BamNode*);
        SR_SWAP(firstType, secondType, int);
    }

    if (firstType == GOOD_ANCHOR)
    {
        switch(secondType)
        {
            case GOOD_ORPHAN:
                return SR_UNIQUE_ORPHAN;
                break;
            case GOOD_SOFT:
                return SR_UNIQUE_SOFT;
                break;
            case GOOD_POOR:
                return SR_UNIQUE_POOR;
                break;
            case GOOD_ANCHOR:
                return SR_UNIQUE_NORMAL;
                break;
            default:
                return SR_OTHER_ALGN_TYPE;
                break;
        }
    }

    return SR_OTHER_ALGN_TYPE;
}

/* 
SR_Bool SR_IsUniqueOrphanPair(SR_BamNode** ppAnchor, SR_BamNode** ppOrphan, double scTolerance, unsigned char minMQ)
{
    int anchorStatus = SR_CheckAlignment(&((*ppAnchor)->alignment), scTolerance, minMQ);
    int orphanStatus = SR_CheckAlignment(&((*ppOrphan)->alignment), scTolerance, minMQ);

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

SR_Status SR_LoadUniquOrphanPairs(SR_BamInStream* pBamInStream, unsigned int threadID, double scTolerance, unsigned char minMQ)
{
    SR_BamNode* pAnchor = NULL;
    SR_BamNode* pOrphan = NULL;

    SR_Status readerStatus = SR_OK;
    SR_Status bufferStatus = SR_OK;
    while ((readerStatus = SR_BamInStreamLoadPair(&pAnchor, &pOrphan, pBamInStream)) == SR_OK)
    {
        if (SR_IsUniqueOrphanPair(&pAnchor, &pOrphan, scTolerance, minMQ))
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
*/

SR_Bool SR_ReadPairFilter(SR_BamNode* pBamNode, const void* filterData)
{
    if ((pBamNode->alignment.core.flag & BAM_FPAIRED) == 0
        || strcmp(bam1_qname(&(pBamNode->alignment)), "*") == 0
        || (pBamNode->alignment.core.flag & SR_READ_PAIR_FMASK) != 0
        || (pBamNode->alignment.core.isize == 0))
    {
        return TRUE;
    }

    // get the statistics of the read pair
    SR_BamPairStats pairStats;
    SR_Status status = SR_LoadPairStats(&pairStats, pBamNode);
    if (status == SR_ERR)
        return TRUE;

    // this is the fragment length distribution
    const SR_FragLenDstrb* pDstrb = (const SR_FragLenDstrb*) filterData;

    // any reads do not have valid read group name will be filtered out
    int32_t readGrpIndex = 0;
    status = SR_FragLenDstrbGetRGIndex(&readGrpIndex, pDstrb, pairStats.RG);
    if (status == SR_ERR)
        return TRUE;
    
    // any reads aligned to different chromosome will be kept as SV candidates
    if (pBamNode->alignment.core.tid != pBamNode->alignment.core.mtid)
    {
        return FALSE;
    }

    // any reads aligned with improper pair mode (orientation) will be kept as SV candidates
    if (!SR_IsValidPairMode(pDstrb, pairStats.pairMode))
    {
        return FALSE;
    }

    // any reads with fragment length at the edge of the fragment length distribution will be kept as SV candidates
    uint32_t lowerCutoffIndex = pDstrb->pHists[readGrpIndex].cutoff[DSTRB_LOWER_CUTOFF];
    uint32_t upperCutoffIndex = pDstrb->pHists[readGrpIndex].cutoff[DSTRB_UPPER_CUTOFF];

    uint32_t lowerCutoff = pDstrb->pHists[readGrpIndex].fragLen[lowerCutoffIndex];
    uint32_t upperCutoff = pDstrb->pHists[readGrpIndex].fragLen[upperCutoffIndex];

    if (pairStats.fragLen < lowerCutoff || pairStats.fragLen > upperCutoff)
        return FALSE;

    // at last, those reads with valid pair mode and proper fragment length will be filtered out
    return TRUE;
}

SR_Status SR_LoadAlgnPairs(SR_BamInStream* pBamInStream, 
                           SR_FragLenDstrb* pDstrb, 
			   bamFile* bam_writer_complete_bam,  // store non-candidate alignments
			   unsigned int threadID, 
			   double scTolerance, // the allowed clips of anchor alignments
			   double maxMismatchRate, // the allowed mismatches of anchor alignments 
			   unsigned char minMQ)
{
    SR_BamNode* pAlgnOne = NULL;
    SR_BamNode* pAlgnTwo = NULL;

    SR_Status readerStatus = SR_OK;
    SR_Status bufferStatus = SR_OK;
    // SR_BamInStreamLoadPair is in utilities/bam/SR_BamInStream.c
    while ((readerStatus = SR_BamInStreamLoadPair(&pAlgnOne, &pAlgnTwo, pBamInStream, bam_writer_complete_bam)) == SR_OK)
    {
        // Load a pair of alignments and store them in pAlgnOne and pAlgnTwo;
	// the position of pAlgnOne is smaller than the one of pAlgnTwo
        SR_AlgnType algnType = SR_GetAlignmentType(&pAlgnOne, &pAlgnTwo, scTolerance, maxMismatchRate, minMQ);
        if ((algnType == SR_UNIQUE_ORPHAN || algnType == SR_UNIQUE_SOFT || algnType == SR_UNIQUE_POOR) && pBamInStream->numThreads > 0)
        {
            #ifdef VERBOSE_DEBUG
	      if (algnType == SR_UNIQUE_ORPHAN) fprintf(stderr, "\n\nSR_UNIQUE_ORPHAN\n");
	      else if (algnType == SR_UNIQUE_SOFT) fprintf(stderr, "\n\nSR_UNIQUE_SOFT\n");
	      else if(algnType == SR_UNIQUE_POOR) fprintf(stderr, "\n\nSR_UNIQUE_POOR\n");
	      else;
	    #endif
	    bufferStatus = SR_BamInStreamPush(pBamInStream, pAlgnOne, threadID);
            bufferStatus = SR_BamInStreamPush(pBamInStream, pAlgnTwo, threadID);

            SR_BamInStreamSetAlgnType(pBamInStream, threadID, algnType);
        }
        else // the pair is not a candidate pair
        {
            if (algnType == SR_UNIQUE_NORMAL && pDstrb != NULL)
            {
                SR_BamPairStats pairStats;
                SR_Status modeStatus = SR_LoadPairStats(&pairStats, pAlgnOne);
                if (modeStatus == SR_OK)
                    SR_FragLenDstrbUpdate(pDstrb, &pairStats);
            }

            // Store alignments in the complete bam
	    if (bam_writer_complete_bam != NULL) {
	      bam_write1(*bam_writer_complete_bam, &(pAlgnOne->alignment));
	      bam_write1(*bam_writer_complete_bam, &(pAlgnTwo->alignment));
	    }

	    SR_BamInStreamRecycle(pBamInStream, pAlgnOne);
            SR_BamInStreamRecycle(pBamInStream, pAlgnTwo);
        }

        if (bufferStatus == SR_FULL)
            break;
    }

    return readerStatus;
}

