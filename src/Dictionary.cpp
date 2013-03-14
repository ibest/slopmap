#include "Dictionary.h"


/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildLibDictionary(string filename) {
    
    std::fstream infile;
    /*Open given file:*/
    infile.open(filename.c_str());
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
    PutLibKmer(contaminant_read, lib_id);
    
    contaminant_read.clear();
    
    return 0;
}


void PutLibKmer(string str, string lib_id)
{
  
  int line_len = str.length();
  /*k_mer_struct - the structure wich holds information about k-mers. This structure is saved in the dictionary.*/
  k_mer_struct k_struct;
  vector<k_mer_struct> rec_id_set;
        
  for(int w=0; w< line_len-KMER_SIZE; w++) {
    
     string seq0;
   
     if((w+KMER_SIZE) < line_len) 
     {
       seq0 = str.substr(w,KMER_SIZE);
     } else 
     {
       seq0 = str.substr(w, line_len-w);
     }
    
     /*Making complement : */
     string seq_complement = MakeSeqComplement(seq0);
                
                
     /*Sequence (read) ID:*/
     k_struct.seq_id = lib_id;
     /*Calculating a position in the read:*/
     k_struct.pos = w;
      
      
     map<string, vector<k_mer_struct> >::iterator it_LibDict = LibDict.find(seq0);
    
     if(it_LibDict == LibDict.end()) 
     {
         rec_id_set.push_back(k_struct);
         LibDict.insert(std::pair<string, vector<k_mer_struct> >(seq0, rec_id_set));
         rec_id_set.clear();
     } else 
     {
         (*it_LibDict).second.push_back(k_struct);
     }
         
         
     it_LibDict = LibDict.find(seq_complement);
    
     if(it_LibDict == LibDict.end()) 
     {
        rec_id_set.push_back(k_struct);
        LibDict.insert(std::pair<string, vector<k_mer_struct> >(seq_complement, rec_id_set));
        rec_id_set.clear();
     } else 
     {
        (*it_LibDict).second.push_back(k_struct);
     }
         
     /*Making reverse complement*/
     reverse( seq_complement.begin(), seq_complement.end() );
     it_LibDict = LibDict.find(seq_complement);
    
     if(it_LibDict == LibDict.end()) 
     {
        rec_id_set.push_back(k_struct);
        LibDict.insert(std::pair<string, vector<k_mer_struct> >(seq_complement, rec_id_set));
        rec_id_set.clear();
     } else 
     {
        (*it_LibDict).second.push_back(k_struct);
     }
         
    }
    
  
}
