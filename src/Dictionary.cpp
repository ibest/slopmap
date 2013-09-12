#include "Dictionary.h"


struct twoBit *twoBitFromDnaSeq(struct dnaSeq *seq, string lib_id)
/* Convert dnaSeq representation in memory to twoBit representation.*/
{
        int ubyteSize = packedSize(seq->size);
        UBYTE *pt;
        struct twoBit *tB;
        DNA last[KMER_SIZE];	/* Holds few bases. */
        DNA *dna;
        int i, end;

        /* Allocate structure and fill in name. */
        tB = (twoBit*)malloc(sizeof(*tB));
        tB->data = (UBYTE*)malloc(sizeof(*(tB->data)) * ubyteSize);
        pt = tB->data;
        //twoBit->name = seq->name;
        //twoBit->size = seq->size;

        /* Convert to 4-bases per byte representation. */
        dna = seq->dna;
        end = seq->size - KMER_SIZE;
        
        /*Calculating a position in the read:*/
        
     
        for (i=0; i<end; i += KMER_SIZE)
        {
                k_mer_struct k_struct;
                k_struct.seq_id = lib_id;
                k_struct.pos = i;
                
                dense_hash_map<UBYTE, vector<k_mer_struct> >::iterator it_LibDict = LibDict.find(packDna(dna+i));
                if(it_LibDict == LibDict.end()) 
                {
                    vector<k_mer_struct> rec_id_set;
                    rec_id_set.push_back(k_struct);
                    LibDict.insert(std::pair<UBYTE, vector<k_mer_struct> >(packDna(dna+i), rec_id_set));
                } else 
                {
                    (*it_LibDict).second.push_back(k_struct);
                }
        }
                
        return tB;
}

/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildLibDictionary2(char* filename) {
    initNtVal();
    UBYTE z = 18446744073709551615;
    LibDict.set_empty_key(z);
    LibDictId.set_empty_key("");
    struct dnaSeq seq;
    ZeroVar(&seq);
    
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
               seq.dna = (char*)contaminant_read.c_str(); seq.size = contaminant_read.length(); seq.name = "Test1";
               twoBitFromDnaSeq(&seq,lib_id);
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
    //cout << contaminant_read << endl;
    /*The last one : */
    seq.dna = (char*)contaminant_read.c_str(); seq.size = contaminant_read.length(); seq.name = "Test1";
    twoBitFromDnaSeq(&seq,lib_id);
    //LibDictId[lib_id] = contaminant_read.length();
    //PutLibKmer(contaminant_read, lib_id);
    
    contaminant_read.clear();
    
    return 0;
}


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
    //cout << contaminant_read << endl;
    /*The last one : */
    LibDictId[lib_id] = contaminant_read.length();
    PutLibKmer(contaminant_read, lib_id);
    
    contaminant_read.clear();
    
    return 0;
}


void PutLibKmer(string str, string lib_id)
{
  /*
  int line_len = str.length();
  //k_mer_struct - the structure wich holds information about k-mers. This structure is saved in the dictionary.
  k_mer_struct k_struct;
  
        
  for(long w=0; w< line_len-KMER_SIZE;  w++) {
      
     string seq0;
   
     if((w+KMER_SIZE) < line_len) 
     {
       seq0 = str.substr(w,KMER_SIZE);
     } else 
     {
       seq0 = str.substr(w, line_len-w);
     }
    //cout << seq0 << endl;
     //Sequence (read) ID:
     k_struct.seq_id = lib_id;
     //Calculating a position in the read:
     k_struct.pos = w;
     //cout << seq0 << endl;
     uint32_t twobit_kmer = dna_number(seq0); 
     //cout << twobit_kmer << endl;
     dense_hash_map<uint32_t, vector<k_mer_struct> >::iterator it_LibDict = LibDict.find(twobit_kmer);
    
     if(it_LibDict == LibDict.end()) 
     {
         vector<k_mer_struct> rec_id_set;
         rec_id_set.push_back(k_struct);
         
         LibDict.insert(std::pair<uint32_t, vector<k_mer_struct> >(twobit_kmer, rec_id_set));
         //LibDict[seq0] = rec_id_set;
     } else 
     {
         (*it_LibDict).second.push_back(k_struct);
     }
         
     //Making reverse complement
     string seq_rev_complement = MakeRevComplement(seq0);
     uint32_t twobit_kmer_revc = dna_number(seq_rev_complement);
     it_LibDict = LibDict.find(twobit_kmer_revc);
    
     if(it_LibDict == LibDict.end()) 
     {
        vector<k_mer_struct> rec_id_set;
        
        rec_id_set.push_back(k_struct);
        //LibDict.insert(std::pair<string, vector<k_mer_struct> >(seq_rev_complement, rec_id_set));
        LibDict[twobit_kmer_revc] = rec_id_set;
     } else 
     {
        (*it_LibDict).second.push_back(k_struct);
     }
         
    }
    
  */
}
