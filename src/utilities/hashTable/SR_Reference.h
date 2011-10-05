/*
 * =====================================================================================
 *
 *       Filename:  SR_Reference.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/03/2011 08:05:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_REFERENCE_H
#define  SR_REFERENCE_H

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>

#include "utilities/common/SR_Types.h"


//===============================
// Type and constant definition
//===============================

// reset the reference object for next reading
#define SR_ReferenceReset(pRef)               \
    do                                        \
    {                                         \
        (pRef)->seqLen = 0;                   \
                                              \
    }while(0) 


// reference header strcture
typedef struct SR_RefHeader
{
    char** names;             // an array contains the name of each chromosome

    void* dict;               // a hash table of reference name and reference ID

    char* md5s;               // an array contains the md5 string of each chromosome

    int64_t* refFilePos;      // an array contains the file offset poisition of each chromosomes in the reference file

    int64_t* htFilePos;       // an array contains the file offset poisition of each chromosomes in the hash table file

    uint32_t numRefs;         // total number of chromosomes 

}SR_RefHeader;

// an object hold the reference sequence of a chromosome
typedef struct SR_Reference
{
    char* sequence;               // reference sequence

    int32_t  id;                  // id of the chromosome

    uint32_t seqLen;              // length of the chromosome

    uint32_t seqCap;              // capacity of reference sequence

}SR_Reference;

// get the reference name given the reference ID
#define SR_RefHeaderGetName(pRefHeader, refID) ((pRefHeader)->names[(refID)])

// get the md5 string (not null terminated) given the reference ID
#define SR_RefHeaderGetMD5(pRefHeader, refID) ((pRefHeader)->md5s + (refID) * MD5_STR_LEN)

//===============================
// Constructors and Destructors
//===============================

// create a new reference object
SR_Reference* SR_ReferenceAlloc(void);

// free an existing reference object
void SR_ReferenceFree(SR_Reference* pRef);

SR_RefHeader* SR_RefHeaderAlloc(void);

void SR_RefHeaderFree(SR_RefHeader* pRefHeader);


//==========================================
// Interface functions related with input
//==========================================

//====================================================================
// function:
//      read the reference header from the reference file
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. refInput: a file pointer to the input reference file
// 
// return:
//      the file offset of the reference header structure in the 
//      reference file. It is used to check the compatibility between
//      the reference file and the hash table file
//====================================================================
int64_t SR_RefHeaderRead(SR_RefHeader* pRefHeader, FILE* refInput);

//===================================================================
// function:
//      get the reference ID given the reference name
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. refName: a reference name
// 
// return:
//      if the reference name is found, return the corresponding
//      reference ID. Otherwise, return -1
//===================================================================
int32_t SR_RefHeaderGetRefID(const SR_RefHeader* pRefHeader, const char* refName);

//====================================================================
// function:
//      jump to a certain chromosome given the reference ID
//
// args:
//      1. refInput: a file pointer to the input reference file
//      2. pRefHeader: a pointer to the reference header structure
//      3. refID: the ID of the reference we want to jump to
// 
// return:
//      SR_OK: successfully jumped
//      SR_ERR: jump failed
//====================================================================
SR_Status SR_ReferenceJump(FILE* refInput, const SR_RefHeader* pRefHeader, int32_t refID);

//====================================================================
// function:
//      read the reference sequence from the input reference file 
//
// args:
//      1. pRef: a pointer to the reference sequence structure
//      2. refInput: a file pointer to the input reference file
// 
//====================================================================
void SR_ReferenceRead(SR_Reference* pRef, FILE* refInput);


//==========================================
// Interface functions related with output
//==========================================

//===================================================================
// function:
//      read the reference sequence in the fasta file line by line, 
//      one chromosome at each time
//
// args:
//      1. pRef: a pointer to the reference structure
//      2. pRefHeader: a pointer to the reference header structure
//      3. faInput: a file pointer to the input fasta file
// 
// return:
//      SR_OK: successfully load the chromosome sequence
//      SR_EOF: reach the end of the fasta file
//      SR_ERR: find an error during loading
//===================================================================
SR_Status SR_ReferenceLoad(SR_Reference* pRef, SR_RefHeader* pRefHeader, FILE* faInput);

//===================================================================
// function:
//      skip the reference sequence with unknown chromosome
//
// args:
//      1. pRefHeader: a pointer to the reference header structure
//      2. faInput: a file pointer to the input fasta file
// 
// return:
//      SR_OK: successfully skipped a chromosome sequence
//      SR_EOF: reach the end of the fasta file
//      SR_ERR: find an error during skipping
//===================================================================
SR_Status SR_ReferenceSkip(SR_RefHeader* pRefHeader, FILE* faInput);

//===================================================================
// function:
//      leave enough space at the beginning of the output reference
//      output file to store the reference header position
//
// args:
//      1. refOutput: a file pointer to the output reference file
//===================================================================
void SR_ReferenceLeaveStart(FILE* refOutput);

//===================================================================
// function:
//      set the reference header position
//
// args:
//      1. refHeaderPos: the reference header position
//      2. refOutput: a file pointer to the output reference file
//===================================================================
void SR_ReferenceSetStart(int64_t refHeaderPos, FILE* refOutput);

//===================================================================
// function:
//      write a reference sequence into the reference output file
//
// args:
//      1. pRef: a pointer to the reference sequence structure
//      2. refOutput: a file pointer to the output reference file
//
// return:
//      the file offset of the current reference sequence
//===================================================================
int64_t SR_ReferenceWrite(const SR_Reference* pRef, FILE* refOutput);

//====================================================================
// function:
//      write the reference header into the output reference file
//
// args:
//      1. pRefHeader: a pointer to the reference sequence structure
//      2. refOutput: a file pointer to the output reference file
//
// return:
//      the file offset of the reference header
//====================================================================
int64_t SR_RefHeaderWrite(const SR_RefHeader* pRefHeader, FILE* refOutput);


#endif  /*SR_REFERENCE_H*/
