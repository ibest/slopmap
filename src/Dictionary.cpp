#include "Dictionary.h"


/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildLibDictionary(char* filename) {
    
    ifstream infile;
    /*Open given file:*/
    if (!infile.is_open())
        infile.open( filename, ifstream::in );
    
    std::string str;
    cout << "Parsing library file " << filename << endl;
    
    /*Length of the read:*/
    int line_len;
    
    /*Counts the number of lines (reads) in parsed file:*/
    long line_cnt = 0;
    
    /*Record (read) ID*/
    string lib_id = "";
    string contaminant_read;
    /*Loop thru all lines in input file:*/
    while ( getline(infile, str) ) {
        
        line_len = str.length();
        if (line_len == 0) {
            continue;
        }
        
        /*Determining a record (read) ID:*/
        remove( str.begin(), str.end(), ' ' );
        if(str.substr(0,1)== ">") {
            if( contaminant_read.length() > 0) 
            {
               LibDictId.insert(std::pair<string, int>(lib_id, contaminant_read.length()));
               PutLibKmer(contaminant_read, lib_id);
               line_cnt = 0;
               contaminant_read.clear();
            }
            lib_id = str.substr(1,str.length()-1);
            cout << lib_id << endl;
            continue;
        }
        
        /*To uppercase:*/
        stoupper(str);
        
        /*Getting the whole read : */
        contaminant_read+=str;
        
        line_cnt+=line_len;
        
    }
    
    infile.close();
    
    /*The last one : */
    LibDictId[lib_id] = contaminant_read.length();
    PutLibKmer(contaminant_read, lib_id);
    
    contaminant_read.clear();
    
    return 0;
}


void PutLibKmer(string str, string lib_id)
{
  
  int line_len = str.length();
  /*k_mer_struct - the structure wich holds information about k-mers. This structure is saved in the dictionary.*/
  k_mer_struct k_struct;
  
        
  for(long w=0; w< line_len-KMER_SIZE; /*w+=DISTANCE*/ ++w) {
    
     string seq0;
   
     if((w+KMER_SIZE) < line_len) 
     {
       seq0 = str.substr(w,KMER_SIZE);
     } else 
     {
       seq0 = str.substr(w, line_len-w);
     }
    
                
                
     /*Sequence (read) ID:*/
     k_struct.seq_id = lib_id;
     /*Calculating a position in the read:*/
     k_struct.pos = w;
     
      
     dense_hash_map<string, vector<k_mer_struct> >::iterator it_LibDict = LibDict.find(seq0);
    
     if(it_LibDict == LibDict.end()) 
     {
         vector<k_mer_struct> rec_id_set;
         rec_id_set.push_back(k_struct);
         //LibDict.insert(std::pair<string, vector<k_mer_struct> >(seq0, rec_id_set));
         LibDict[seq0] = rec_id_set;
     } else 
     {
         (*it_LibDict).second.push_back(k_struct);
     }
         
     /*Making reverse complement*/
     string seq_rev_complement = MakeRevComplement(seq0);
     if (seq0 == seq_rev_complement ) 
     {
         cout << seq0 << " " << seq_rev_complement << endl;
     }
     it_LibDict = LibDict.find(seq_rev_complement);
    
     if(it_LibDict == LibDict.end()) 
     {
        vector<k_mer_struct> rec_id_set;
        rec_id_set.push_back(k_struct);
        //LibDict.insert(std::pair<string, vector<k_mer_struct> >(seq_rev_complement, rec_id_set));
        LibDict[seq_rev_complement] = rec_id_set;
     } else 
     {
        (*it_LibDict).second.push_back(k_struct);
     }
         
    }
    
  
}
