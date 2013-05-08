#include "KMerRoutine.h"


LibHitData CheckForLib(string seq)
//unsigned int CheckForLib(string seq)
{
    //Sampling with some special frequency:
    string read = string(seq);
    int line_len = seq.length();
    stoupper(read);
    //cout <<  line_len << endl;
    //Number of k-mers in the read:
    unsigned short knum = ceil((double)(line_len - KMER_SIZE + 1)/DISTANCE);
    unsigned short numhits = 0;
    
    LibHitData most_similar;
    most_similar.start_pos = -1;
    int i = 0;
    
    while( i+KMER_SIZE< line_len ) 
    {
       string kmer_str;
       
       kmer_str = read.substr(i,KMER_SIZE);
       
       it_LibDict = LibDict.find(kmer_str);
       if(it_LibDict != LibDict.end()) 
       {
           /*Got hit:*/
           if((*it_LibDict).second.size() == 1) 
           {
               numhits += 1;
               
               if(((double)numhits/(double)knum) >= similarity_threshold) 
               {
                   most_similar.lib_id = (*it_LibDict).second[0].seq_id;
                   most_similar.start_pos = (*it_LibDict).second[0].pos;
                   most_similar.end_pos = (*it_LibDict).second[0].pos;
                   break;
               }
           }
       }
       
       i = i+DISTANCE;   
   }
    
   return most_similar;
}

vector<LibHitData> Similarity(vector <HitData> matches, short read_len) 
{
    map<string, int* > scores; //Represent a pair: <LibId,counterb>
    map<string, int*>::iterator it_scores;
    
    for(unsigned short i=0; i<matches.size(); ++i)
    {
        for(unsigned short j=0; j<matches[i].kmers.size(); ++j) 
        {
            it_scores = scores.find(matches[i].kmers[j].seq_id);
            if(it_scores == scores.end())
            {
                int *a = NULL;
                a = new int[3];
                a[0] = 1;
                a[1] = matches[i].kmers[j].pos;
                a[2] = matches[i].kmers[j].pos;
                scores.insert(std::pair<string,int* >(matches[i].kmers[j].seq_id,a)); //LibId, start pos, end pos
            }
            else
            {
                int tmp = (*it_scores).second[1];
                (*it_scores).second[0] = (*it_scores).second[0] + 1;
                (*it_scores).second[1] = tmp;
                (*it_scores).second[2] = matches[i].kmers[j].pos;
            }
            
        }
        //cout << endl;
    }
    
    //Now we have to find a maximum:
    //Build a second map with keys that are the first's map values and values as the first map's keys.
    map<int, string > tmp_scores;
    it_scores = scores.begin();
    while(it_scores != scores.end())
    {
        tmp_scores.insert(std::pair<int,string>((*it_scores).second[0], (*it_scores).first));
        it_scores++;
    }
    
    vector<LibHitData> most_similar;
    
    //Now threshold:
    if(tmp_scores.size() > 0)
    {
        map<int, string>::reverse_iterator it_tmp_scores;
        it_tmp_scores = tmp_scores.rbegin();
        double score = 0.0;
        while(it_tmp_scores != tmp_scores.rend())
        {
                score = ((double)(*it_tmp_scores).first/(double)(floor((read_len-KMER_SIZE)/DISTANCE)+1))*100;
                if(score >= similarity_threshold)
                {
                    
                    LibHitData most_similar_lib;
                    //cout << (double)(*it_tmp_scores).first << (*it_tmp_scores).second << endl;
                    most_similar_lib.lib_id = (*it_tmp_scores).second;
                    //cout << (*scores.find(most_similar_lib.lib_id)).second[1] << endl;
                    most_similar_lib.start_pos = (*scores.find(most_similar_lib.lib_id)).second[1];
                    most_similar_lib.end_pos = (*scores.find(most_similar_lib.lib_id)).second[2];
                    most_similar.push_back(most_similar_lib);
                }
                else
                {
                    break;
                }
                it_tmp_scores++;
                
        }
    }
    
    
    tmp_scores.clear();
    scores.clear();
    
    return most_similar;
}

