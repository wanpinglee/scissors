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
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "utilities/md5.h"
#include "hasher/common/SR_Error.h"
#include "SR_Reference.h"

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
static unsigned char ProcessHeaderLine(const char* buff)
{
    // we only accept two formats of header
    // 1) the chromosome ID should closely followed by the '>' character.
    // 2) the chromosome ID should closely followed by the ">chr" string.
    // the end of chromosome ID is either detected with a space character, a tab character, a new line character or a null character.

    const char* header = buff + 1;
    if (strncmp(header, "chr", 3) == 0)
        header = buff + 4;

    // the ID of a chromosome should be less than 4 characters
    char chrBuff[4];

    for(unsigned int i = 0; i != 3; ++i)
    {
        if (header[i] == ' ' || header[i] == '\t' || header[i] == '\n' || header[i] == '\0')
        {
            chrBuff[i] = '\0';
            break;
        }

        chrBuff[i] = header[i];
    }

    chrBuff[3] = '\0';

    // by default chrX is denoted as 23
    // chrY is denoted as 24
    // chrMT is denoted as 25
    // the rest chromosome will use a arabic number as its ID
    // 0 IS NOT ALLOWED TO BE USED AS A CHROMOSOME ID. 0 IS USED AS AN ERROR FLAG.

    if (strcmp(chrBuff, "X") == 0)
        return X;
    else if (strcmp(chrBuff, "Y") == 0)
        return Y;
    else if (strcmp(chrBuff, "MT") == 0)
        return MT;
    else
        return (unsigned char) atoi(chrBuff);
}

// calculate the md5 sum value of a reference sequence
// and store the result in a reference object
static void SetMd5(SR_Reference* reference)
{
    unsigned char MD5[MD5_CHECKSUM_LEN];

    memset(MD5, 0, MD5_CHECKSUM_LEN);

    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char*) reference->sequence, reference->length);
    MD5Final(MD5, &context);

    char* md5String = reference->md5;
    for (unsigned int i = 0; i != MD5_CHECKSUM_LEN; ++i)
    {
        sprintf(md5String, "%02x", MD5[i]);
        md5String += 2;
    }
}


// create a new reference object
SR_Reference* SR_ReferenceAlloc(uint32_t capacity)
{
    SR_Reference* newRef = (SR_Reference*) malloc(sizeof(SR_Reference));
    if (newRef == NULL)
        SR_ErrSys("ERROR: Not enough memory for a reference object.\n");

    if (capacity == 0)
    {
        SR_ErrMsg("WARNING: Capacity of reference sequence should be greater than zero. A default value %d will be used.\n", DEFAULT_REF_CAPACITY);
        capacity = DEFAULT_REF_CAPACITY;
    }

    newRef->sequence = (char*) malloc(sizeof(char) * capacity);
    if (newRef->sequence == NULL)
        SR_ErrSys("ERROR: Not enough memory for the storage of sequence in a reference object.\n");

    newRef->chr = 0;
    newRef->length = 0;
    newRef->capacity = capacity;
    newRef->md5[MD5_STR_LEN] = '\0';

    return newRef;
}

// free an existing reference object
void SR_ReferenceFree(SR_Reference* reference)
{
    if (reference != NULL)
    {
        free(reference->sequence);

        free(reference);
    }
}

// read the reference sequence in the fasta file line by line, one chromosome at each time
SR_Bool SR_ReferenceLoad(SR_Reference* reference, unsigned char* nextChr, FILE* faInput)
{
    char buff[MAX_REF_LINE];

    while (fgets(buff, MAX_REF_LINE, faInput) != NULL && buff[0] != '>')
    {
        // actual length of the reference line
        // excluding the starting and trailing space and the new line character
        unsigned short length = 0;
        const char* refLine = ProcessRefLine(&length, buff);

        // we skip the following processes if we get an empty line
        if (refLine == NULL)
            continue;

        // we shouldn't get this step if we are handling human reference
        // the default capacity(300Mbp) should be enough even for the largest chromosome in human genome
        if (length + reference->length > reference->capacity)
        {
            reference->capacity *= 2;
            reference->sequence = (char*) realloc(reference->sequence, sizeof(char) * reference->capacity);
            if (reference->sequence == NULL) 
                SR_ErrQuit("ERROR: Not enough memory for the storage of sequence in the reference object.\n");
        }

        // copy the current line into the reference object
        char* cpAddr= reference->sequence + reference->length;
        strncpy(cpAddr, refLine, length);
        reference->length += length;
    }

    if (reference->length > 0)
        SetMd5(reference);

    if (buff[0] == '>')
        *nextChr = ProcessHeaderLine(buff);
    else if (feof(faInput))
        return FALSE;
    else
        SR_ErrSys("ERROR: Found an error when reading the fasta file.");

    return TRUE;
}

// skip the reference sequence with unknown chromosome ID
SR_Bool SR_ReferenceSkip(unsigned char* nextChr, FILE* faInput)
{

    char buff[MAX_REF_LINE];

    while (fgets(buff, MAX_REF_LINE, faInput) != NULL && buff[0] != '>')
        continue;

    if (buff[0] == '>')
        *nextChr = ProcessHeaderLine(buff);
    else if (feof(faInput))
        return FALSE;
    else
        SR_ErrSys("ERROR: Found an error when reading the fasta file.");

    return TRUE;
}

// write the reference sequence into a output file in the binary format
off_t SR_ReferenceWrite(FILE* refOutput, const SR_Reference* reference)
{
    off_t fileOffset = ftello(refOutput);
    if (fileOffset == -1)
        SR_ErrSys("ERROR: Cannot get the offset of current file.\n");

    size_t writeSize = 0;

    writeSize = fwrite(&(reference->chr), sizeof(unsigned char), 1, refOutput);
    if (writeSize != 1)
        SR_ErrSys("ERROR: Cannot write chromosome ID into the reference file.\n");

    writeSize = fwrite(reference->md5, sizeof(char), MD5_STR_LEN, refOutput);
    if (writeSize != MD5_STR_LEN)
        SR_ErrSys("ERROR: Cannot write md5 string into the reference file.\n");

    writeSize = fwrite(&(reference->length), sizeof(uint32_t), 1, refOutput);
    if (writeSize != 1)
        SR_ErrSys("ERROR: Cannot write chromosome length into the reference file.\n");

    writeSize = fwrite(reference->sequence, sizeof(char), reference->length, refOutput);
    if (writeSize != reference->length)
        SR_ErrSys("ERROR: Cannot write chromosome sequence into the reference file.\n");

    fflush(refOutput);

    return fileOffset;
}

// read the reference sequence from an input file in the binary format
SR_Bool SR_ReferenceRead(SR_Reference* reference, FILE* refInput)
{
    size_t readSize = 0;

    readSize = fread(&(reference->chr), sizeof(unsigned char), 1, refInput);
    if (readSize != 1)
    {
        if (feof(refInput))
            return FALSE;
        else
            SR_ErrQuit("ERROR: Cannot read the reference input file due to an error.\n");
    }

    readSize = fread(reference->md5, sizeof(char), MD5_STR_LEN, refInput);
    if (readSize != MD5_STR_LEN)
        SR_ErrSys("ERROR: Cannot read md5 string from the reference file.\n");

    readSize = fread(&(reference->length), sizeof(uint32_t), 1, refInput);
    if (readSize != 1)
        SR_ErrSys("ERROR: Cannot read chromosome length from the reference file.\n");

    if (reference->length > reference->capacity)
    {
        reference->capacity = reference->length;

        free(reference->sequence);
        reference->sequence = (char*) malloc(sizeof(char) * reference->capacity);
        if (reference->sequence == NULL)
            SR_ErrQuit("ERROR: Not enough memory for reference sequence.\n");
    }

    readSize = fread(reference->sequence, sizeof(char), reference->length, refInput);
    if (readSize != reference->length)
        SR_ErrSys("ERROR: Cannot read chromosome sequence from the reference file.\n");

    return TRUE;
}
