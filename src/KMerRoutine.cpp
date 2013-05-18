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
       it_LibDict = LibDict.find(kmer_str);
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



void CheckForLib2(string seq) 
{
    /*
    string read = string(seq);
    int line_len = seq.length();
    stoupper(read);
    int ii = 0;
    comp_TODO = ceil((double)(line_len - KMER_SIZE));
    while( ii+KMER_SIZE< line_len )  // fill inQueue will rubbish data since this isn't an actual computation...
    {
        string kmer_str;
        kmer_str = read.substr(ii,KMER_SIZE);
        inQueue.push_back(kmer_str);
        ii++;
   }
   
    
   
    for (unsigned long i=0; i<16; i++) // start the threads
    {
        pthread_t *tId( new pthread_t );   threadIdList.push_back(tId);
        thread_detail Y; Y.num=i; thread_table.push_back(Y);
        int rc( pthread_create( tId, NULL, workerThread, (void *)(&(thread_table.back() )) ) );
        if (rc) { std::cout<<"ERROR; return code from pthread_create() "<<comp_START<<"\n"; std::cout.flush();
              exit(-1); }
   }
  
   // now we wait for the threads to terminate, perhaps updating the screen with info as we go. 
   
   while (comp_DONE != comp_TODO)
   {
      // poll the queue to get a status update on computation
      pthread_mutex_lock(&outQueueLock);
      comp_DONE = outQueue.size();
      pthread_mutex_unlock(&outQueueLock);
   } // big while loop

   // call join to kill all worker threads
   std::list< pthread_t* >::iterator i(threadIdList.begin());
   while (i!=threadIdList.end())
   {
    if (pthread_join( *(*i), NULL)!=0) { std::cout<<"Thread join error!\n"; exit(1); }
    delete (*i);
    threadIdList.erase(i++);  
   }
   //std::cout<<"\n";
   
   
    // clean-up
   inQueue.clear();
   outQueue.clear();     
     */   
    // pthread_exit(NULL);
}
