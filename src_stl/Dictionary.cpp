#include "Dictionary.h"


struct twoBit *twoBitFromDnaSeq(struct dnaSeq *seq, string lib_id)
/* Convert dnaSeq representation in memory to twoBit representation.*/
{
        int ubyteSize = packedSize(seq->size);
//        UBYTE *pt;
        struct twoBit *tB;
 //       DNA last[KMER_SIZE];	/* Holds few bases. */
        DNA *dna;
        int i, end;

        /* Allocate structure and fill in name. */
        tB = (twoBit*)malloc(sizeof(*tB));
        tB->data = (UBYTE*)malloc(sizeof(*(tB->data)) * ubyteSize);
//        pt = tB->data;
        //twoBit->name = seq->name;
        //twoBit->size = seq->size;

        /* Convert to 4-bases per byte representation. */
        dna = seq->dna;
        end = seq->size - KMER_SIZE;
        
        /*Calculating a position in the read:*/
        map<UBYTE, vector<k_mer_struct> > LibDict;
        
     
        for (i=0; i<end; ++i)//KMER_SIZE)
        {
                k_mer_struct k_struct;
                k_struct.seq_id = lib_id;
                k_struct.pos = i;
                
                //dense_hash_map<UBYTE, vector<k_mer_struct> >::iterator it_LibDict = LibDict.find(packDna(dna+i));
                map<UBYTE, vector<k_mer_struct> >::iterator it_LibDict = LibDict.find(packDna(dna+i));
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
        
        dict_holder.push_back(LibDict);
        
        return tB;
}

/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildLibDictionary2(char* filename) {
    initNtVal();
    
    
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
    
    if ( (rep_file = fopen((char*)rep_file_name.c_str(), "w")) == NULL ) 
    {
        fprintf(stderr,"[err] Could not open file '%s' for reading.\n", (char*)rep_file_name.c_str());
        exit(1);
    }
    
    
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
               //LibDictId.insert(std::pair<string, int>(lib_id, 0));
               seq.dna = (char*)contaminant_read.c_str(); seq.size = contaminant_read.length();
               
               twoBitFromDnaSeq(&seq,lib_id);
               line_cnt = 0;
               contaminant_read.clear();
               //LibDict.clear();
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
    //LibDictId.insert(std::pair<string, int>(lib_id, 0));
    seq.dna = (char*)contaminant_read.c_str(); seq.size = contaminant_read.length();
    twoBitFromDnaSeq(&seq,lib_id);
    
    if(illumina_pe_flag)
    {
        IlluminaDynamic();
    } else if(roche_flag)
    {
       Roche454Dynamic();
    } else if(illumina_se_flag)
    {
      IlluminaDynamicSE();
    }
    
    contaminant_read.clear();
    //LibDict.clear();
    
    pe_output_file1.close();
    pe_output_file2.close();
    
    cout << "====================Done====================\n";  
    fclose(rep_file);
    
    return 0;
}