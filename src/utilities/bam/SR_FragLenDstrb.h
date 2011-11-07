/*
 * =====================================================================================
 *
 *       Filename:  SR_FragLenDstrb.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/07/2011 03:06:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  SR_FRAGLENDSTRB_H
#define  SR_FRAGLENDSTRB_H

#include "outsources/samtools/bam.h"
#include "utilities/common/SR_Types.h"
#include "SR_BamHeader.h"
#include "SR_BamMemPool.h"
// #include "SR_BamPairAux.h"

//===============================
// Type and constant definition
//===============================

enum SR_FragLenDstrbCutoffIndex
{
    DSTRB_LOWER_CUTOFF = 0,

    DSTRB_UPPER_CUTOFF = 1
};

// a map used to map the pair mode into its corresponding number
// negative value means invalid mode
static const int SR_PairModeMap[16] = { 
                                          -1, -1, 0, 1,

                                          -1, -1, 2, 3,

                                          4, 5, -1, -1,

                                          6, 7, -1, -1
                                      };


static const int SR_PairModeSetMap[64] =  { 
                                               0, 0, 0, 1, 0, 0, 0, 1,
                                               0, 0, 1, 0, 0, 1, 0, 0,
                                               0, 1, 0, 0, 0, 0, 1, 0,
                                               1, 0, 0, 0, 1, 0, 0, 0,
                                               0, 0, 0, 1, 0, 0, 0, 1,
                                               0, 1, 0, 0, 0, 0, 1, 0,
                                               0, 0, 1, 0, 0, 1, 0, 0,
                                               1, 0, 0, 0, 1, 0, 0, 0
                                          };


// the object used to hold the basic statistics of a pair of alignments
typedef struct SR_BamPairStats
{
    const char* RG;           // name of the read group

    unsigned int fragLen;     // fragment length of the pair

    uint16_t pairMode;        // orientation mode of the pair

}SR_BamPairStats;

// the object used to hold the fragment length histogram of a given read group
typedef struct SR_FragLenHist
{
    void* rawHist[NUM_ALLOWED_HIST];                  // raw fragment length histogram. this is a hash table for histogram building

    uint32_t* fragLen[NUM_ALLOWED_HIST];              // array of the fragment length

    double* cdf[NUM_ALLOWED_HIST];                    // array of the cumulative probability at a given fragment length

    double mean[NUM_ALLOWED_HIST];                    // mean of the histogram

    double median[NUM_ALLOWED_HIST];                  // median of the histogram

    double stdev[NUM_ALLOWED_HIST];                   // standard deviation of the histogram

    uint32_t size[NUM_ALLOWED_HIST];                  // number of unique fragment length

    uint32_t cutoff[NUM_ALLOWED_HIST][2];             // cutoff of the probability

    uint64_t modeCount[NUM_ALLOWED_HIST + 1];         // total counts of a histogram

}SR_FragLenHist;

// the object used to hold the fragment length distribution of all the read groups in a bam file
typedef struct SR_FragLenDstrb
{
    char** pReadGrpNames;                          // read group name array

    void* pReadGrpHash;                            // hash used to map the read group name to read group index

    SR_FragLenHist* pHists;                        // array of fragment length histograms

    int8_t validModeMap[NUM_TOTAL_PAIR_MODE];      // map the pair mode to its corresponding histogram, invalid pair mode will get a negative value
    
    int8_t validMode[NUM_ALLOWED_PAIR_MODE];       // the valid pair modes
    
    uint8_t numPairMode;                           // number of valid pair modes given by the user

    uint32_t size;                                 // number of read groups found in the bam file

    uint32_t capacity;                             // capacity of the read group name array

    SR_Bool hasRG;                                 // boolean variable to indicate if we find any read group names or not in the current bam file

}SR_FragLenDstrb;


//===============================
// Constructors and Destructors
//===============================

SR_FragLenDstrb* SR_FragLenDstrbAlloc(uint32_t capacity);

void SR_FragLenDstrbFree(SR_FragLenDstrb* pDstrb);


//======================
// Inline functions
//======================

//==================================================================
// function:
//      get the pair orientation mode from a bam alignment
//
// args:
//      1. pBamNode: a pointer to a bam node structure
//
// return:
//      the mode of the pair
//==================================================================
SR_PairMode SR_GetPairMode(const SR_BamNode* pBamNode);

//=======================================================================
// function:
//      get the basic statistics from a pair of alignments
//
// args:
//      1. pPairStats: a poiter to the pair statistics object
//      2. ppUpAlgn: a pointer to a bam node object for the alignment 
//                   with smaller coordinate
//      3. ppDownAlgn: a pointer to a bam node object for the alignment 
//                     with greater coordinate
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_LoadPairStats(SR_BamPairStats* pPairStats, const SR_BamNode* pBamNode);


void SR_FragLenDstrbSetPairMode(SR_FragLenDstrb* pDstrb, const int8_t* pValidPairMode, uint8_t numPairMode);

//=======================================================================
// function:
//      retrieve the read group information from the bam header and
//      initialize the fragment length distribution with it
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. pBamHeader: a pointer to a bam header object
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_FragLenDstrbSetRG(SR_FragLenDstrb* pDstrb, const SR_BamHeader* pBamHeader);

//=======================================================================
// function:
//      given a read group name, find its corresponding read group ID
//
// args:
//      1. pReadGrpIndex: a pointer to the read group index variable
//                        this is the return value
//      2. pDstrb: a pointer to a fragment length distribution object
//      3. pReadGrpName: a pointer to the read group name
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_FragLenDstrbGetRGIndex(int32_t* pReadGrpIndex, const SR_FragLenDstrb* pDstrb, const char* pRreadGrpName);

//=======================================================================
// function:
//      update the fragment length distribution with the information
//      from a new read pair
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. pPairStats: a pointer to a read pair statistic object
//
// return:
//      if an error occurs, return SR_ERR; else return SR_OK
//========================================================================
SR_Status SR_FragLenDstrbUpdate(SR_FragLenDstrb* pDstrb, const SR_BamPairStats* pPairStats);

//=======================================================================
// function:
//      after all the read pairs are processed, raw histogram will be
//      transferred into mature one. its mean, median and standard
//      deviation will be calculated.
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//========================================================================
void SR_FragLenDstrbFinalize(SR_FragLenDstrb* pDstrb);

//=======================================================================
// function:
//      set the probability cutoff for all the histograms
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. cutoff: the fragment length probability cutoff
//========================================================================
void SR_FragLenDstrbSetCutoff(SR_FragLenDstrb* pDstrb, double cutoff);

//=======================================================================
// function:
//      write the fragment length distribution into a file
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. dstrbOutput: output stream
//========================================================================
void SR_FragLenDstrbWrite(const SR_FragLenDstrb* pDstrb, FILE* dstrbOutput);

//=======================================================================
// function:
//      read the fragment length distribution from a file
//
// args:
//      1. pDstrb: a pointer to a fragment length distribution object
//      2. dstrbOutput: input stream
//========================================================================
void SR_FragLenDstrbRead(SR_FragLenDstrb* pDstrb, FILE* dstrbInput);

#endif  /*SR_FRAGLENDSTRB_H*/
