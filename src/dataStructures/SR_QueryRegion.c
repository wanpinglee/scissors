/*
 * =====================================================================================
 *
 *       Filename:  SR_QueryRegion.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/06/2011 12:42:45 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>

#include "hashTable/common/SR_Error.h"
#include "hashTable/common/SR_Utilities.h"
#include "SR_QueryRegion.h"

// map used to transfer the 4-bit representation of a nucleotide into the ascii representation
static const char SR_BASE_MAP[16] = {'N', 'A', 'C', 'N', 'G','N','N','N','T','N','N','N','N','N','N','N'};

SR_QueryRegion* SR_QueryRegionAlloc(void)
{
    SR_QueryRegion* pNewRegion = (SR_QueryRegion*) malloc(sizeof(SR_QueryRegion));
    if (pNewRegion == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a orphanSeq region object");

    // be caureful, the actual data are not stored in these two pointers
    // they are stored in the bam aux data structure
    pNewRegion->pAnchor = NULL;
    pNewRegion->pOrphan = NULL;

    pNewRegion->orphanSeq = NULL;
    pNewRegion->isOrphanInversed = FALSE;
    pNewRegion->capacity = 0;

    pNewRegion->closeRefBegin = 0;
    pNewRegion->closeRefEnd = 0;
    pNewRegion->farRefBegin = 0;
    pNewRegion->farRefEnd = 0;

    return pNewRegion;
}

void SR_QueryRegionFree(SR_QueryRegion* pQueryRegion)
{
    if (pQueryRegion != NULL)
    {
        free(pQueryRegion->orphanSeq);
        free(pQueryRegion);
    }
}


void SR_QueryRegionLoadSeq(SR_QueryRegion* pQueryRegion)
{
    // if we don't have enough space for the orphan sequence, we need to expand current storage space
    if (pQueryRegion->capacity < pQueryRegion->pOrphan->core.l_qseq)
    {
        pQueryRegion->capacity = pQueryRegion->orphanSeq != NULL ? pQueryRegion->pOrphan->core.l_qseq * 2 : pQueryRegion->pOrphan->core.l_qseq;
        free(pQueryRegion->orphanSeq);

        pQueryRegion->orphanSeq = (char*) malloc(sizeof(char) * pQueryRegion->capacity);
        if (pQueryRegion->orphanSeq == NULL)
            SR_ErrQuit("ERROR: Not enough memory for a orphanSeq string.");
    }

    uint8_t* seq = bam1_seq(pQueryRegion->pOrphan);
    for (unsigned int i = 0; i != pQueryRegion->pOrphan->core.l_qseq; ++i)
    {
        pQueryRegion->orphanSeq[i] = SR_BASE_MAP[bam1_seqi(seq, i)];
    }
}

void SR_QueryRegionSetSeq(SR_QueryRegion* pQueryRegion, SR_SeqAction action)
{
    if (action == SR_INVERSE || action == SR_REVERSE_COMP)
    {
        for (unsigned int i = 0, j = pQueryRegion->pOrphan->core.l_qseq - 1; i < j; ++i, --j)
        {
	    SR_SWAP(pQueryRegion->orphanSeq[i], pQueryRegion->orphanSeq[j], char);
        }
    }

    if (action == SR_COMP || action == SR_REVERSE_COMP)
    {
        for (unsigned int i = 0; i != pQueryRegion->pOrphan->core.l_qseq; ++i)
        {
            switch(pQueryRegion->orphanSeq[i])
            {
                case 'A':
                    pQueryRegion->orphanSeq[i] = 'T';
                    break;
                case 'C':
                    pQueryRegion->orphanSeq[i] = 'G';
                    break;
                case 'G':
                    pQueryRegion->orphanSeq[i] = 'C';
                    break;
                case 'T':
                    pQueryRegion->orphanSeq[i] = 'A';
                    break;
                default:
                    pQueryRegion->orphanSeq[i] = 'N';
                    break;
            }
        }
    }
}

SR_Bool SR_QueryRegionSetRange(SR_QueryRegion* pQueryRegion, const SR_SearchArgs* pSearchArgs, uint32_t refLen, SR_Direction direction)
{
    if (direction == SR_DOWNSTREAM) // the position of the search region is greater than that of the anchor mate
    {
	// calculate the rightmost coordinate of an alignment on the reference genome
	uint32_t onRefEnd = bam_calend(&(pQueryRegion->pAnchor->core), bam1_cigar(pQueryRegion->pAnchor));

        pQueryRegion->closeRefBegin = onRefEnd + 1 + pSearchArgs->fragLen - pSearchArgs->closeRange / 2;
        // out of range (close/far begin is larger than the reference length of current chr)
        if (pQueryRegion->closeRefBegin >= refLen)
            return FALSE;

        // in this case, the far reference begins equals to the close reference begin
        pQueryRegion->farRefBegin = pQueryRegion->closeRefBegin;

        pQueryRegion->closeRefEnd = pQueryRegion->closeRefBegin + pSearchArgs->closeRange - 1;
        if (pQueryRegion->closeRefEnd >= refLen)
            pQueryRegion->closeRefEnd = refLen - 1;

        // search region is too small to hold a read
        if (pQueryRegion->closeRefEnd < pQueryRegion->closeRefBegin + pQueryRegion->pAnchor->core.l_qseq - 1)
            return FALSE;

        pQueryRegion->farRefEnd = pQueryRegion->farRefBegin + pSearchArgs->farRange - 1;
        if (pQueryRegion->farRefEnd >= refLen)
            pQueryRegion->farRefEnd = refLen - 1;
    }
    else if (direction == SR_UPSTREAM) // the position of the search region is less than that of the anchor mate
    {
        if (pQueryRegion->pAnchor->core.pos  + pSearchArgs->closeRange / 2 > 1 + pSearchArgs->fragLen)
            pQueryRegion->closeRefEnd = pQueryRegion->pAnchor->core.pos + pSearchArgs->closeRange / 2 - 1 - pSearchArgs->fragLen;
        else
            return FALSE; // out of range (close/far end is smaller than 0)


        // search region is too small to hold a read
        if (pQueryRegion->closeRefEnd - 0 + 1 < pQueryRegion->pAnchor->core.l_qseq)
            return FALSE;

        // in this case, the far reference end equals to the close reference end
        pQueryRegion->farRefEnd = pQueryRegion->closeRefEnd;

        if (pQueryRegion->closeRefEnd + 1 >= pSearchArgs->closeRange)
            pQueryRegion->closeRefBegin = pQueryRegion->closeRefEnd + 1 - pSearchArgs->closeRange;
        else
            pQueryRegion->closeRefBegin = 0;

        if (pQueryRegion->farRefEnd + 1 >= pSearchArgs->farRange)
            pQueryRegion->farRefBegin = pQueryRegion->farRefEnd + 1 - pSearchArgs->farRange;
        else
            pQueryRegion->farRefBegin = 0;
    }

    return TRUE;
}
