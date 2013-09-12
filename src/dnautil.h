/* 
 * File:   dnautil.h
 * Author: ilya
 *
 * Created on September 11, 2013, 4:01 PM
 */

#ifndef DNAUTIL_H
#define	DNAUTIL_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <iostream>
#include <fstream>


/* Numerical values for bases. */
#define MASKED_BASE_BIT 8
#define T_BASE_VAL 0
#define U_BASE_VAL 0
#define C_BASE_VAL 1
#define A_BASE_VAL 2
#define G_BASE_VAL 3
#define N_BASE_VAL 4   /* Used in 1/2 byte representation. */
/* Some other type synonyms */
//#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
#define UBYTE unsigned long long
#define BYTE signed char      /* Wants to be signed 8 bits. */
#define UWORD unsigned short  /* Wants to be unsigned 16 bits. */
#define WORD short	      /* Wants to be signed 16 bits. */
#define bits64 unsigned long long  /* Wants to be unsigned 64 bits. */
#define bits32 unsigned       /* Wants to be unsigned 32 bits. */
#define bits16 unsigned short /* Wants to be unsigned 16 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
#define signed32 int	      /* Wants to be signed 32 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */

#define boolean int

/* How big is this array? */
#define ArraySize(a) (sizeof(a)/sizeof((a)[0]))

#ifndef min
#define min(a,b) ( (a) < (b) ? (a) : (b) )
/* Return min of a and b. */
#endif

/* inline functions: To declare a function inline, place the entire function
 * in a header file and prefix it with the INLINE macro.  If used with a
 * compiler that doesn't support inline, change the INLINE marco to be simply
 * `static'.
 */
#ifndef INLINE
#define INLINE static inline
#endif

#define NEEDMEM_LIMIT 500000000

#define AllocVar(pt) (pt = needMem(sizeof(*pt)))
/* Shortcut to allocating a single variable on the heap and
 * assigning pointer to it. */

#define AllocArray(pt, size) (pt = needLargeZeroedMem(sizeof(*pt) * (size)))

INLINE void zeroBytes(void *vpt, int count)
/* fill a specified area of memory with zeroes */
{
  memset(vpt, '\0', count);
}

#define ZeroVar(v) zeroBytes(v, sizeof(*v))

typedef char DNA;


extern short KMER_SIZE;

struct twoBit
/* Two bit representation of DNA. */
{
    struct twoBit *next;	/* Next sequence in list */
    char *name;			/* Name of sequence. */
    UBYTE *data;		/* DNA at two bits per base. */
    bits32 size;		/* Size of this sequence. */
    bits32 nBlockCount;		/* Count of blocks of Ns. */
    bits32 *nStarts;		/* Starts of blocks of Ns. */
    bits32 *nSizes;		/* Sizes of blocks of Ns. */
    bits32 maskBlockCount;	/* Count of masked blocks. */
    bits32 *maskStarts;		/* Starts of masked regions. */
    bits32 *maskSizes;		/* Sizes of masked regions. */
    bits32 reserved;		/* Reserved for future expansion. */
};

struct dnaSeq
/* A dna sequence in one-character per base format. */
{
    struct dnaSeq *next;  /* Next in list. */
    char *name;           /* Name of sequence. */
    DNA *dna;             /* Sequence base by base. */
    int size;             /* Size of sequence. */
    //Bits* mask;           /* Repeat mask (optional) */
};

//static void initNtVal();
void initNtVal();
int packedSize(int unpackedSize);
//static int packedSize(int unpackedSize);
UBYTE packDna4(DNA *in);
UBYTE packDna(DNA *in);
struct twoBit *twoBitFromDnaSeq(struct dnaSeq *seq);

#endif	/* DNAUTIL_H */

