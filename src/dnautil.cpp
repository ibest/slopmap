#include "dnautil.h"

int ntVal[256];
int ntValLower[256];	/* NT values only for lower case. */
int ntValUpper[256];	/* NT values only for upper case. */
int ntVal5[256];
int ntValNoN[256]; /* Like ntVal, but with T_BASE_VAL in place of -1 for nonexistent ones. */
DNA valToNt[(N_BASE_VAL|MASKED_BASE_BIT)+1];

/* convert tables for bit-4 indicating masked */
int ntValMasked[256];
DNA valToNtMasked[256];
boolean inittedNtVal = false;


//static void initNtVal()
void initNtVal()
{
        if (!inittedNtVal)
        {
                unsigned int i;
                for (i=0; i<ArraySize(ntVal); i++)
                {
                        ntValUpper[i] = ntValLower[i] = ntVal[i] = -1;
                        ntValNoN[i] = T_BASE_VAL;
                        if (isspace(i) || isdigit(i))
                                ntVal5[i] = ntValMasked[i] = -1;
                        else
                        {
                                ntVal5[i] = N_BASE_VAL;
                                ntValMasked[i] = (islower(i) ? (N_BASE_VAL|MASKED_BASE_BIT) : N_BASE_VAL);
                        }
                }
                ntVal5['t'] = ntVal5['T'] = ntValNoN['t'] = ntValNoN['T'] = ntVal['t'] = ntVal['T'] = 
                ntValLower['t'] = ntValUpper['T'] = T_BASE_VAL;
                ntVal5['u'] = ntVal5['U'] = ntValNoN['u'] = ntValNoN['U'] = ntVal['u'] = ntVal['U'] = 
                ntValLower['u'] = ntValUpper['U'] = U_BASE_VAL;
                ntVal5['c'] = ntVal5['C'] = ntValNoN['c'] = ntValNoN['C'] = ntVal['c'] = ntVal['C'] = 
                ntValLower['c'] = ntValUpper['C'] = C_BASE_VAL;
                ntVal5['a'] = ntVal5['A'] = ntValNoN['a'] = ntValNoN['A'] = ntVal['a'] = ntVal['A'] = 
                ntValLower['a'] = ntValUpper['A'] = A_BASE_VAL;
                ntVal5['g'] = ntVal5['G'] = ntValNoN['g'] = ntValNoN['G'] = ntVal['g'] = ntVal['G'] = 
                ntValLower['g'] = ntValUpper['G'] = G_BASE_VAL;

                valToNt[T_BASE_VAL] = valToNt[T_BASE_VAL|MASKED_BASE_BIT] = 't';
                valToNt[C_BASE_VAL] = valToNt[C_BASE_VAL|MASKED_BASE_BIT] = 'c';
                valToNt[A_BASE_VAL] = valToNt[A_BASE_VAL|MASKED_BASE_BIT] = 'a';
                valToNt[G_BASE_VAL] = valToNt[G_BASE_VAL|MASKED_BASE_BIT] = 'g';
                valToNt[N_BASE_VAL] = valToNt[N_BASE_VAL|MASKED_BASE_BIT] = 'n';


                inittedNtVal = true;
        }
}

//static int packedSize(int unpackedSize)
int packedSize(int unpackedSize)
/* Return size when packed, rounding up. */
{
        return ((unpackedSize + 3) >> 2);
}


UBYTE packDna4(DNA *in)
/* Pack 4 bases into a UBYTE */
{
  UBYTE out = 0;
  int count = 4;
  int bVal;
  while (--count >= 0) {
    bVal = ntValNoN[(int)*in++];
    out <<= 2;
    out += bVal;
  }
  return out;
}

UBYTE packDna(DNA *in)
/* Pack KMER_SIZE bases into a UBYTE */
{
  UBYTE out = 0;
  int count = KMER_SIZE;
  int bVal;
  while (--count >= 0) {
    bVal = ntValNoN[(int)*in++];
    out <<= 2;
    out += bVal;
    //printf("%c",*in);
  }
  //printf("\n");
  //printf("%d\n",out);
  return out;
}

