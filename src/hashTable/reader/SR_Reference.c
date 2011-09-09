/*
 * =====================================================================================
 *
 *       Filename:  SR_SR_Reference.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/03/2011 08:06:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "utilities/md5.h"
#include "samtools/khash.h"
#include "hashTable/common/SR_Error.h"
#include "SR_Reference.h"


//===============================
// Type and constant definition
//===============================

// maximum number of character will be load in a line from the fasta file
#define MAX_REF_LINE 1024

// number of characters can be held in a reference object
// this value assure that the largest chromosome in the human reference, chromsome 1, can be
// load into the object without any reallocation.
#define DEFAULT_REF_CAP 300000000

// the default start chromosome ID
#define DEFAULT_START_CHR 1

// the default number of chromosomes
#define DEFAULT_NUM_CHR 100

// initialize the hash table for reference name
KHASH_MAP_INIT_STR(refName, int32_t);


//===================
// Static methods
//===================

// process a line of reference sequence
static const char* ProcessRefLine(unsigned short* len, char* buff)
{
    char* iter = buff;
    const char* refLine = buff;
    *len = 0;

    while (*refLine == ' ' || *refLine == '\t')
    {
        if (*refLine == '\n' || *refLine == '\0')
        {
            refLine = NULL;
            return refLine;
        }

        ++iter;
        ++refLine;
    }

    while (*iter != '\n' && *iter != ' ' && *iter != '\t' && *iter != '\0')
    {
        *iter = toupper(*iter);
        ++iter;
        ++(*len);
    }

    return refLine;
}

// process the header line in the fasta file to get the ID for the next chromosome
static void SR_RefHeaderSetName(SR_RefHeader* pRefHeader, const char* buff)
{
    // we only accept two formats of header
    // 1) the chromosome ID should closely followed by the '>' character.
    // 2) the chromosome ID should closely followed by the ">chr" string.
    // the end of chromosome ID is either detected with a space character, a tab character, a new line character or a null character.


    // skip the '>' character
    const char* header = buff + 1;
    while (*header != ' ' && *header != '\t' 
           && *header != '\n' && *header != '\0')
    {
        ++header;
    }

    unsigned int headerLen = header - buff - 1;
    if (headerLen != 0)
    {
        pRefHeader->names[pRefHeader->numRefs] = (char*) malloc((headerLen + 1) * sizeof(char));
        if (pRefHeader->names[pRefHeader->numRefs] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the reference name.\n");

        pRefHeader->names[pRefHeader->numRefs][headerLen] = '\0';
        strncpy(pRefHeader->names[pRefHeader->numRefs], buff + 1, headerLen);
    }
    else
        pRefHeader->names[pRefHeader->numRefs] = NULL;
}

static void SR_RefHeaderSetMd5(SR_RefHeader* pRefHeader, SR_Reference* pRef)
{
    pRef->id = pRefHeader->numRefs;

    unsigned char MD5[MD5_CHECKSUM_LEN];
    memset(MD5, 0, MD5_CHECKSUM_LEN);

    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char*) pRef->sequence, pRef->seqLen);
    MD5Final(MD5, &context);

    char* md5String = pRefHeader->md5s + MD5_STR_LEN * pRefHeader->numRefs;
    for (unsigned int i = 0; i != MD5_CHECKSUM_LEN; ++i)
    {
        sprintf(md5String, "%02X", MD5[i]);
        md5String += 2;
    }

    ++(pRefHeader->numRefs);
}


//===============================
// Constructors and Destructors
//===============================

// create a new reference object
SR_Reference* SR_ReferenceAlloc(void)
{
    SR_Reference* newRef = (SR_Reference*) malloc(sizeof(SR_Reference));
    if (newRef == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a reference object.\n");

    newRef->sequence = (char*) malloc(sizeof(char) * DEFAULT_REF_CAP);
    if (newRef->sequence == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of sequence in a reference object.\n");

    newRef->id = 0;
    newRef->seqLen = 0;
    newRef->seqCap = DEFAULT_REF_CAP;

    return newRef;
}

// free an existing reference object
void SR_ReferenceFree(SR_Reference* pRef)
{
    if (pRef != NULL)
    {
        free(pRef->sequence);
        free(pRef);
    }
}

SR_RefHeader* SR_RefHeaderAlloc(void)
{
    SR_RefHeader* pRefHeader = (SR_RefHeader*) calloc(1, sizeof(SR_RefHeader));
    if (pRefHeader == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a reference header object.\n");

    pRefHeader->names = (char**) calloc(DEFAULT_NUM_CHR, sizeof(char*));
    if (pRefHeader->names == NULL)
        SR_ErrQuit("ERROR: Not enough memory for reference names in a reference header object.\n");

    pRefHeader->dict = (void*) kh_init(refName);

    pRefHeader->md5s = (char*) calloc(DEFAULT_NUM_CHR * MD5_STR_LEN, sizeof(char));
    if (pRefHeader->md5s == NULL)
        SR_ErrQuit("ERROR: Not enough memory for MD5 strings in a reference header object.\n");

    pRefHeader->refFilePos = (int64_t*) calloc(DEFAULT_NUM_CHR, sizeof(int64_t));
    if (pRefHeader->refFilePos == NULL)
        SR_ErrQuit("ERROR: Not enough memory for reference file positions in a reference header object.\n");

    pRefHeader->htFilePos = (int64_t*) calloc(DEFAULT_NUM_CHR, sizeof(int64_t));
    if (pRefHeader->htFilePos == NULL)
        SR_ErrQuit("ERROR: Not enough memory for hash table file positions in a reference header object.\n");

    pRefHeader->numRefs = 0;

    return pRefHeader;
}

void SR_RefHeaderFree(SR_RefHeader* pRefHeader)
{
    if (pRefHeader != NULL)
    {
        kh_destroy(refName, pRefHeader->dict);

        if (pRefHeader->names != NULL)
        {
            for (unsigned int i = 0; i != pRefHeader->numRefs; ++i)
                free(pRefHeader->names[i]);

            free(pRefHeader->names);
        }

        free(pRefHeader->md5s);
        free(pRefHeader->refFilePos);
        free(pRefHeader->htFilePos);

        free(pRefHeader);
    }
}


//==========================================
// Interface functions related with input
//==========================================

// read the reference header from the reference file
int64_t SR_RefHeaderRead(SR_RefHeader* pRefHeader, FILE* refInput)
{
    size_t readSize = 0;
    int64_t refHeaderPos = 0;
    int64_t refStart = 0;

    readSize = fread(&refHeaderPos, sizeof(int64_t), 1, refInput);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the reference header position from the reference file.\n");

    if ((refStart = ftello(refInput)) < 0)
        SR_ErrQuit("ERROR: Cannot get the offset of current file.\n");

    if (fseeko(refInput, refHeaderPos, SEEK_SET) != 0)
        SR_ErrQuit("ERROR: Cannot seek in the reference file.\n");

    readSize = fread(&(pRefHeader->numRefs), sizeof(uint32_t), 1, refInput);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the number of chromosomes from the reference file.\n");

    int khRet = 0;
    khiter_t iter;
    khash_t(refName)* hash = pRefHeader->dict;
    for (unsigned int i = 0; i != pRefHeader->numRefs; ++i)
    {
        uint32_t nameLen = 0;
        readSize = fread(&nameLen, sizeof(uint32_t), 1, refInput);
        if (readSize != 1)
            SR_ErrQuit("ERROR: Cannot read the reference name length from the reference file.\n");

        pRefHeader->names[i] = (char*) malloc(sizeof(char) * (nameLen + 1));
        if (pRefHeader->names[i] == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the reference name.\n");

        pRefHeader->names[i][nameLen] = '\0';
        readSize = fread(pRefHeader->names[i], sizeof(char), nameLen, refInput);
        if (readSize != nameLen)
            SR_ErrQuit("ERROR: Cannot read the reference name from the reference file.\n");

        iter = kh_put(refName, hash, pRefHeader->names[i], &khRet);
        kh_value(hash, iter) = i;
    }

    readSize = fread(pRefHeader->md5s, sizeof(char), MD5_STR_LEN * pRefHeader->numRefs, refInput);
    if (readSize != MD5_STR_LEN * pRefHeader->numRefs)
        SR_ErrQuit("ERROR: Cannot read the md5 strings from the reference file.\n");

    readSize = fread(pRefHeader->refFilePos, sizeof(int64_t), pRefHeader->numRefs, refInput);
    if (readSize !=  pRefHeader->numRefs)
        SR_ErrQuit("ERROR: Cannot read the offset of chromosomes from the reference file.\n");

    readSize = fread(pRefHeader->htFilePos, sizeof(int64_t), pRefHeader->numRefs, refInput);
    if (readSize != pRefHeader->numRefs)
        SR_ErrQuit("ERROR: Cannot read the offset of hash table from the reference file.\n");

    if (fseeko(refInput, refStart, SEEK_SET) != 0)
        SR_ErrQuit("ERROR: Cannot seek in the reference file.\n");

    return refHeaderPos;
}

// get the reference ID given the reference name
int32_t SR_RefHeaderGetRefID(const SR_RefHeader* pRefHeader, const char* refName)
{
    khiter_t iter;
    khash_t(refName)* hash = (khash_t(refName)*) pRefHeader->dict;
    iter = kh_get(refName, hash, refName);

    return iter == kh_end(hash)? -1 : kh_value(hash, iter);
}

// jump to a certain chromosome given the reference ID
SR_Status SR_ReferenceJump(FILE* refInput, const SR_RefHeader* pRefHeader, int32_t refID)
{
    int64_t jumpPos = pRefHeader->refFilePos[refID];

    if (fseeko(refInput, jumpPos, SEEK_SET) != 0)
        return SR_ERR;

    return SR_OK;
}

// read the reference sequence from the input reference file 
void SR_ReferenceRead(SR_Reference* pRef, FILE* refInput)
{
    size_t readSize = 0;

    readSize = fread(&(pRef->id), sizeof(int32_t), 1, refInput);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read the pRef input file due to an error.\n");

    readSize = fread(&(pRef->seqLen), sizeof(uint32_t), 1, refInput);
    if (readSize != 1)
        SR_ErrQuit("ERROR: Cannot read chromosome length from the reference file.\n");

    if (pRef->seqLen > pRef->seqCap)
    {
        pRef->seqCap = pRef->seqLen;

        free(pRef->sequence);
        pRef->sequence = (char*) malloc(sizeof(char) * pRef->seqCap);
        if (pRef->sequence == NULL)
            SR_ErrQuit("ERROR: Not enough memory for reference sequence.\n");
    }

    readSize = fread(pRef->sequence, sizeof(char), pRef->seqLen, refInput);
    if (readSize != pRef->seqLen)
        SR_ErrQuit("ERROR: Cannot read chromosome sequence from the reference file.\n");
}


//==========================================
// Interface functions related with output
//==========================================

// read the reference sequence in the fasta file line by line, one chromosome at each time
SR_Status SR_ReferenceLoad(SR_Reference* pRef, SR_RefHeader* pRefHeader, FILE* faInput)
{
    char buff[MAX_REF_LINE];

    while (fgets(buff, MAX_REF_LINE, faInput) != NULL && buff[0] != '>')
    {
        // actual length of the pRef line
        // excluding the starting and trailing space and the new line character
        unsigned short length = 0;
        const char* refLine = ProcessRefLine(&length, buff);

        // we skip the following processes if we get an empty line
        if (refLine == NULL)
            continue;

        // we shouldn't get this step if we are handling human pRef
        // the default capacity(300Mbp) should be enough even for the largest chromosome in human genome
        if (length + pRef->seqLen > pRef->seqCap)
        {
            pRef->seqCap *= 2;
            pRef->sequence = (char*) realloc(pRef->sequence, sizeof(char) * pRef->seqCap);
            if (pRef->sequence == NULL) 
                SR_ErrQuit("ERROR: Not enough memory for the storage of sequence in the pRef object.\n");
        }

        // copy the current line into the pRef object
        char* cpAddr= pRef->sequence + pRef->seqLen;
        strncpy(cpAddr, refLine, length);
        pRef->seqLen += length;
    }

    if (pRef->seqLen > 0)
        SR_RefHeaderSetMd5(pRefHeader, pRef);
    else
    {
        free(pRefHeader->names[pRefHeader->numRefs]);
        pRefHeader->names[pRefHeader->numRefs] = NULL;
    }

    if (buff[0] == '>')
        SR_RefHeaderSetName(pRefHeader, buff);
    else if (feof(faInput))
        return SR_EOF;
    else
        return SR_ERR;

    return SR_OK;
}

// skip the reference sequence with unknown chromosome ID
SR_Status SR_ReferenceSkip(SR_RefHeader* pRefHeader, FILE* faInput)
{
    char buff[MAX_REF_LINE];

    while (fgets(buff, MAX_REF_LINE, faInput) != NULL && buff[0] != '>')
        continue;

    if (buff[0] == '>')
        SR_RefHeaderSetName(pRefHeader, buff);
    else if (feof(faInput))
        return SR_EOF;
    else
        return SR_ERR;

    return SR_OK;
}

// leave enough space at the beginning of the output reference
// output file to store the reference header position
void SR_ReferenceLeaveStart(FILE* refOutput)
{
    int64_t emptyOffset = 0;
    size_t writeSize = 0;

    writeSize = fwrite(&emptyOffset, sizeof(int64_t), 1, refOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the offset of the reference header into the reference file.\n");
}

// set the reference header position
void SR_ReferenceSetStart(int64_t refHeaderPos, FILE* refOutput)
{
    if (fseeko(refOutput, 0, SEEK_SET) != 0)
        SR_ErrSys("ERROR: Cannot seek in the reference output file\n");

    size_t writeSize = 0;

    writeSize = fwrite(&refHeaderPos, sizeof(int64_t), 1, refOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the offset of the reference header into the reference file.\n");
}

// write a reference sequence into the reference output file
int64_t SR_ReferenceWrite(const SR_Reference* pRef, FILE* refOutput)
{
    int64_t offset = ftello(refOutput);
    if (offset < 0)
        SR_ErrQuit("ERROR: Cannot get the offset of current file.\n");

    size_t writeSize = 0;

    writeSize = fwrite(&(pRef->id), sizeof(int32_t), 1, refOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write chromosome ID into the reference file.\n");

    writeSize = fwrite(&(pRef->seqLen), sizeof(uint32_t), 1, refOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write chromosome length into the reference file.\n");

    writeSize = fwrite(pRef->sequence, sizeof(char), pRef->seqLen, refOutput);
    if (writeSize != pRef->seqLen)
        SR_ErrQuit("ERROR: Cannot write chromosome sequence into the reference file.\n");

    fflush(refOutput);

    return offset;
}

// write the reference header into the output reference file
int64_t SR_RefHeaderWrite(const SR_RefHeader* pRefHeader, FILE* refOutput)
{
    int64_t offset = ftello(refOutput);
    if (offset < 0)
        SR_ErrQuit("ERROR: Cannot get the offset of current file.\n");

    size_t writeSize = 0;

    writeSize = fwrite(&(pRefHeader->numRefs), sizeof(uint32_t), 1, refOutput);
    if (writeSize != 1)
        SR_ErrQuit("ERROR: Cannot write the total number of chromosomes into the reference file.\n");

    for (unsigned int i = 0; i != pRefHeader->numRefs; ++i)
    {
        uint32_t nameLen = strlen(pRefHeader->names[i]);
        writeSize = fwrite(&nameLen, sizeof(uint32_t), 1, refOutput);
        if (writeSize != 1)
            SR_ErrQuit("ERROR: Cannot write the length of reference name into the reference file.\n");

        writeSize = fwrite(pRefHeader->names[i], sizeof(char),  nameLen, refOutput);
        if (writeSize != nameLen)
            SR_ErrQuit("ERROR: Cannot write the reference name into the reference file.\n");
    }

    writeSize = fwrite(pRefHeader->md5s, sizeof(char), MD5_STR_LEN * pRefHeader->numRefs, refOutput);
    if (writeSize != MD5_STR_LEN * pRefHeader->numRefs)
        SR_ErrQuit("ERROR: Cannot write the MD5 strings into the reference file.\n");

    writeSize = fwrite(pRefHeader->refFilePos, sizeof(int64_t), pRefHeader->numRefs, refOutput);
    if (writeSize != pRefHeader->numRefs)
        SR_ErrQuit("ERROR: Cannot write the file offsets of the reference into the reference file.\n");

    writeSize = fwrite(pRefHeader->htFilePos, sizeof(int64_t), pRefHeader->numRefs, refOutput);
    if (writeSize != pRefHeader->numRefs)
        SR_ErrQuit("ERROR: Cannot write the file offsets of the hash table into the reference file.\n");

    fflush(refOutput);

    return offset;
}

