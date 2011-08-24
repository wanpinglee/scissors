/*
 * =====================================================================================
 *
 *       Filename:  SR_BamInStream.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/18/2011 05:39:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_BAMINSTREAM_H
#define  SR_BAMINSTREAM_H

#include <stdlib.h>

#include "samtools/bam.h"
#include "hashTable/common/SR_Types.h"


//===============================
// Type and constant definition
//===============================

// structure holding related bam input variables
typedef struct SR_BamInStreamPrvt SR_BamInStream;


//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(const char* bamFilename);

void SR_BamInStreamFree(SR_BamInStream* pBamInStream);


//======================
// Interface functions
//======================

//===============================================================
// function:
//      jump to a certain chromosome in a bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. refID : the reference ID we want to jump to
// 
// return:
//      if jumping succeeds, return SR_OK; if not, return SR_ERR
//=============================================================== 
SR_Status SR_BamInStreamJump(SR_BamInStream* pBamInStream, int32_t refID);

//================================================================
// function:
//      read the header of a bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      if reading succeeds, return a pointer to a header
//      structure; if error is found, return NULL
//
// discussion: 
//      The file position indicator must be placed at the 
//      beginning of the file. Upon success, the position 
//      indicator will be set at the start of the first alignment
//================================================================ 
bam_header_t* SR_BamInStreamReadHeader(SR_BamInStream* pBamInStream);

//================================================================
// function:
//      read an alignment from the bam file
//
// args:
//      1. pAlignment: a pointer to an alignment
//      2. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      if reading succeeds, return SR_OK; if reach the end of
//      file, return SR_EOF; if get error, return SR_ERR
//
// discussion: 
//      The file position indicator must be placed at the 
//      beginning of an alignment. Upon success, the position 
//      indicator will be set at the start of the next alignment
//================================================================ 
SR_Status SR_BamInStreamRead(bam1_t* pAlignment, SR_BamInStream* pBamInStream);

//================================================================
// function:
//      load a unique-orphan pair from the bam file
//
// args:
//      1. ppAnchor: a pointer of pointer to the anchor mate
//      2. ppOrphan: a pointer of pointer to the orphan mate
//      3. pBamInStream : a pointer to an bam instream structure
//
// return:
//      if we get a unique-orphan pair, return SR_OK; if we reach
//      the end of file, return SR_EOF; if we finish the current
//      chromosome, return SR_OUT_OF_RANGE; else, return SR_ERR
//================================================================
SR_Status SR_BamInStreamGetPair(bam1_t** ppAnchor, bam1_t** ppOrphan, SR_BamInStream* pBamReadAux);


#endif  /*SR_BAMINSTREAM_H*/
