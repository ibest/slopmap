#include "KMerRoutine.h"

// the muticies, protectors of the shared resources
extern pthread_mutex_t coutLock;
extern pthread_mutex_t inQueueLock;
extern pthread_mutex_t outQueueLock;
// the shared data
extern std::list< std::string > inQueue;
extern std::list< std::string > outQueue;

 extern unsigned long comp_DONE; 
 extern unsigned long comp_START;
 extern unsigned long comp_TODO;
 
 

LibHitData CheckForLib(string seq)
//unsigned int CheckForLib(string seq)
{
    //Sampling with some special frequency:
    string read = string(seq);
    int line_len = seq.length();
    stoupper(read);
    //cout <<  line_len << endl;
    //Number of k-mers in the read:
    //unsigned short knum = ceil((double)(line_len - KMER_SIZE + 1)/DISTANCE);
//    unsigned short knum = ceil((double)(line_len/KMER_SIZE));
    unsigned short numhits = 0;
    
    LibHitData most_similar;
    most_similar.start_pos = -1;
    int i = 0;
    short i0 = 0;
    short total_len = 0;
    while( i+KMER_SIZE< line_len ) 
    {
       string kmer_str;
       
       kmer_str = read.substr(i,KMER_SIZE);
       //pthread_mutex_lock(&outQueueLock);
       uint32_t kmer_str_bit = dna_number(kmer_str);
       it_LibDict = LibDict.find(kmer_str_bit);
       //pthread_mutex_unlock(&outQueueLock);
       if(it_LibDict != LibDict.end()) 
       {
           if((*it_LibDict).second.size() == 1) 
           {
               numhits += 1;
               
               if(i - i0 >= KMER_SIZE) 
               {
                   total_len += KMER_SIZE;
               }
               else
               {
                   total_len += i - i0;
               }
               
               if(((double)total_len/(double)line_len) >= similarity_threshold) 
               {
                   most_similar.lib_id = (*it_LibDict).second[0].seq_id;
                   most_similar.start_pos = (*it_LibDict).second[0].pos;
                   most_similar.end_pos = (*it_LibDict).second[0].pos;
                   break;
               }
               
               i0 = i;
               
           }
       }
       
       i = i + DISTANCE;   
   }
    
   return most_similar;
}



LibHitData CheckForLib2(string seq)
//unsigned int CheckForLib(string seq)
{
    
    
    //Sampling with some special frequency:
    string read = string(seq);
    int line_len = seq.length();
    stoupper(read);
    //cout <<  line_len << endl;
    //Number of k-mers in the read:
    //unsigned short knum = ceil((double)(line_len - KMER_SIZE + 1)/DISTANCE);
//    unsigned short knum = ceil((double)(line_len/KMER_SIZE));
    unsigned short numhits = 0;
    
    LibHitData most_similar;
    most_similar.start_pos = -1;
    short i0 = 0;
    short total_len = 0;
    
    struct dnaSeq _seq;
    ZeroVar(&_seq);
    _seq.dna = (char*)seq.c_str(); _seq.size = line_len; _seq.name = "Test1";
    
    int ubyteSize = packedSize(_seq.size);
    UBYTE *pt;
    struct twoBit *tB;
    DNA last[KMER_SIZE];	/* Holds few bases. */
    DNA *dna;
    int i, end;
   
    /* Convert to 4-bases per byte representation. */
    dna = _seq.dna;
    end = _seq.size - KMER_SIZE;
    
    /* Allocate structure and fill in name. */
    tB = (twoBit*)malloc(sizeof(*tB));
    tB->data = (UBYTE*)malloc(sizeof(*(tB->data)) * ubyteSize);
    pt = tB->data;
        
    for (i=0; i<end; i += DISTANCE)
    {
        //twoBit->name = seq->name;
        //twoBit->size = seq->size;
        
        UBYTE kmer_str_bit = packDna(dna+i);
        it_LibDict = LibDict.find(kmer_str_bit);
        //pthread_mutex_unlock(&outQueueLock);
        if(it_LibDict != LibDict.end()) 
        {
           if((*it_LibDict).second.size() == 1) 
           {
               numhits += 1;
               
               if(i - i0 >= KMER_SIZE) 
               {
                   total_len += KMER_SIZE;
               }
               else
               {
                   total_len += i - i0;
               }
               
               if(((double)total_len/(double)line_len) >= similarity_threshold) 
               {
                   most_similar.lib_id = (*it_LibDict).second[0].seq_id;
                   most_similar.start_pos = (*it_LibDict).second[0].pos;
                   most_similar.end_pos = (*it_LibDict).second[0].pos;
                   break;
               }
               
               i0 = i;
               
           }
       }
       
   }
    
   return most_similar;
}

