/*
 * =====================================================================================
 *
 *       Filename:  SR_Types.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/08/2011 21:43:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_TYPES_H
#define  SR_TYPES_H

typedef enum
{
    FALSE = 0,
    TRUE  = 1

}SR_Bool;

typedef enum
{
    SR_OK  =           0,
    SR_EOF =          -1,
    SR_ERR =          -2,
    SR_OUT_OF_RANGE = -100

}SR_Status;

typedef enum
{
    SR_A = 1,
    SR_C = 2,
    SR_G = 4,
    SR_T = 8,
    SR_N = 15
}SR_Base;

// search direction for query region relative to the positon of anchor mate
typedef enum
{
    SR_UPSTREAM,     // search the query region uptream to the anchor mate (smaller coordinate)

    SR_DOWNSTREAM    // search the query region downstream to the anchor mate (larger coordinate)

}SR_Direction;

// action applied on a certain DNA sequence
typedef enum
{
    SR_FORWARD,         // leave the sequence as the format read from bam file

    SR_REVERSE_COMP,    // change the sequence to its reverse complement format

    SR_INVERSE          // change the sequence to its inverse format

}SR_SeqAction;

typedef SR_SeqAction SR_Strand;

#define MD5_CHECKSUM_LEN 16

#define MD5_STR_LEN 32

#define X 23

#define Y 24

#define MT 25

#define INVALID_CHR_ID 0

#define MAX_HASH_SIZE 12

#define SR_EMPTY 0

#endif  /*SR_TYPES_H*/
