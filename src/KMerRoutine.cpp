#include "KMerRoutine.h"


string CheckForLib(string seq)
{
    //Sampling with some special frequency:
    string read = string(seq);
    int line_len = seq.length();
    stoupper(read);
    //cout <<  line_len << endl;
    bool stop_flag = false;
    
    vector <HitData> matches;
    int i = 0;
    
    string lib_id = "";
    
    while( i< line_len ) 
    {
       string kmer_str;
       
       kmer_str = read.substr(i,KMER_SIZE);
       
       it_LibDict = LibDict.find(kmer_str);
       if(it_LibDict != LibDict.end()) 
       {
           /*Got hit:*/
           HitData t_match_struct;
           t_match_struct.pos = i;
           t_match_struct.kmers = (*it_LibDict).second;
           t_match_struct.k_mer_string = kmer_str;
           matches.push_back(t_match_struct);
           
           if(!mode_flag) {
                if((int)matches.size() == NUM_CONSEQUITIVE_HITS) 
                {
                        stop_flag = true;
                        lib_id = t_match_struct.kmers[0].seq_id;
                        
                        break;
                }
           }
           else 
           {
               if((int)matches.size() == NUM_HITS) 
                {
                        stop_flag = true;
                        lib_id = t_match_struct.kmers[0].seq_id;
                        break;
                }
           }
       }
               
       if(stop_flag == true) break;   
              
       i = i+ KMER_SIZE+DISTANCE;   
   }
    
   matches.clear();
    
   return lib_id;
}

