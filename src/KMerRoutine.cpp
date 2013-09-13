#include "KMerRoutine.h"


struct LibHitData CheckForLib2(string *seq, dense_hash_map<UBYTE, vector<k_mer_struct> > *LibDict)
{
    //Sampling with some special frequency:
    string read = string(*seq);
    int line_len = seq->length();
    stoupper(read);
    
    unsigned short numhits = 0;
    
    struct LibHitData most_similar;
    //most_similar = (LibHitData*)malloc(sizeof(*most_similar));
    most_similar.start_pos = -1;
    short i0 = 0;
    short total_len = 0;
    
    struct dnaSeq _seq;
    ZeroVar(&_seq);
    _seq.dna = (char*)seq->c_str(); _seq.size = line_len;
    
 //   int ubyteSize = packedSize(_seq.size);
 //   UBYTE *pt;
 //   struct twoBit *tB;
//    DNA last[KMER_SIZE];	/* Holds few bases. */
    DNA *dna;
    int i, end;
   
    /* Convert to 4-bases per byte representation. */
    dna = _seq.dna;
    end = _seq.size - KMER_SIZE;
    
    /* Allocate structure and fill in name. */
    //tB = (twoBit*)malloc(sizeof(*tB));
    //tB->data = (UBYTE*)malloc(sizeof(*(tB->data)) * ubyteSize);
    //pt = tB->data;
        
    for (i=0; i<end; i += DISTANCE)
    {
        UBYTE kmer_str_bit = packDna(dna+i);
        
        dense_hash_map<UBYTE, vector<k_mer_struct> >::iterator it_LibDict;
        it_LibDict = LibDict->find(kmer_str_bit);
        
        if(it_LibDict != LibDict->end()) 
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
           
           i0 = i;
           
           if(((double)total_len/(double)line_len) >= similarity_threshold) 
           {
              most_similar.lib_id = (*it_LibDict).second[0].seq_id;
              most_similar.start_pos = (*it_LibDict).second[0].pos;
              most_similar.end_pos = i;
              break;
           }
       }
       
   }
    
    //free(tB->data);
    //free(tB);
    
   return most_similar;
}

void WritePEFile(fstream &pe_output_file, Read *read)
{
    pe_output_file << read->illumina_readID << endl;
    pe_output_file << read->read << endl;
    pe_output_file << '+' << endl;
    pe_output_file << read->illumina_quality_string << endl;
}

void Roche454Dynamic()
{
    fstream rep_file;
    rep_file.open(rep_file_name.c_str(),ios::out);
    
    for(int i=0; i<(int)roche_names.size(); ++i)
    {
            cout << "Parsing file: " << roche_names[i] << "..." << endl;
            string b;
            //If SFF format is given -> process it
            if( string(roche_names[i]).substr( strlen(roche_names[i])-3, 3 ) == "sff" ) 
            {
                cout << "File is in SFF format.\n" ;
                
                sff_common_header h;
                sff_read_header rh;
                sff_read_data rd;
                FILE *sff_fp;

                if ( (sff_fp = fopen(roche_names[i], "r")) == NULL ) 
                {
                        fprintf(stderr,
                                "[err] Could not open file '%s' for reading.\n", roche_names[i]);
                        exit(1);
                }
               
                read_sff_common_header(sff_fp, &h);
                verify_sff_common_header((char*)PRG_NAME, (char*)SFF_FILE_VERSION, &h);

                int left_clip = 0, right_clip = 0;
                char *bases;
                uint8_t *quality;
                
                int numreads = (int) h.nreads;
                
                for (int i = 0; i < numreads; i++) 
                { 
                    read_sff_read_header(sff_fp, &rh);
                    read_sff_read_data(sff_fp, &rd, h.flow_len, rh.nbases);

                    get_clip_values(rh, 0, &left_clip, &right_clip);
                    
                    // create bases string 
                    bases = get_read_bases(rd, left_clip, right_clip);

                    // create quality array 
                    quality = get_read_quality_values(rd, left_clip, right_clip);
                    b = string(bases);
                    for(int k =0; k<dict_holder.size(); k++) {
                        struct LibHitData match = CheckForLib2(&b,&dict_holder[k]);
                        if(match.start_pos != -1) 
                        {
                                rep_file << match.lib_id << "\t" << match.start_pos << "\t" << match.end_pos << "\t" << rh.name << "\t" << bases << "\t" << quality << "\n";
                        }
                    }           
                    
                    free(bases);
                    free(quality);
                    free_sff_read_header(&rh);
                    free_sff_read_data(&rd);

                }

                fclose(sff_fp);

               
            } 
            else if(string(roche_names[i]).substr( strlen(roche_names[i])-5, 5 ) == "fastq") 
            {
               //FASTQ file given. Process it.
               cout << "File is in FASTQ format, starting conversion...\n" ;
               
               string line, bases, quality, readID;
               int ii = 0;
               igzstream in(roche_names[i]);
               while ( getline(in,line) ) 
               {
                   if(ii==0) readID = line; /*Read ID*/
                
                   if(ii==1) bases = line; /*actual data*/
                
                   if(ii==3) 
                   {
                       for(int k =0; k<dict_holder.size(); k++) {
                                struct LibHitData match = CheckForLib2(&b,&dict_holder[k]);
                                if(match.start_pos != -1) 
                                {
                                        rep_file << match.lib_id << "\t" << match.start_pos << "\t" << match.end_pos << "\t" << readID << "\t" << bases << "\t" << line << "\n";
                                }
                       }
                       
                       ii = 0;
                   }
                   ii++;
               }
               
               in.close();
               
            }
    }
    
    rep_file.close();
}

//Dynamic Illumina: does not need space to store reads:
void IlluminaDynamic()
{
    
    vector<string> record_block1, record_block2;
    
    if(!pe_output_file1.is_open()) {
                                        pe_output_file1.open( pe_output_filename1.c_str(), ios::out );
                                }
                                if(!pe_output_file2.is_open()) {
                                        pe_output_file2.open( pe_output_filename2.c_str(), ios::out );
                                }
                   
    
    for(int jj=0; jj<(int)pe1_names.size(); ++jj)
    {
        int ii = 0;
        
        std::string line1, line2;
        igzstream in1( pe1_names[jj] ); //for R1
        igzstream in2( pe2_names[jj] ); //for R2
        
        cout << "Processing files: " << pe1_names[jj] << ", " << pe2_names[jj] << "\n";
        
        while ( getline(in1,line1) && getline(in2,line2) )
        {
            /*Read ID*/
            if(ii==0) 
            {
                //Check for order
                vector <string> fields1, fields2;
                split_str( line1, fields1, " " );
                split_str( line2, fields2, " " );
                if( (fields1[0] != fields2[0] ) && !old_style_illumina_flag)
                {
                    cout << "Warning: read IDs do not match in input files: PE1-> " << pe1_names[jj] << ", PE2-> " << pe2_names[jj] << endl;
                }
                        
                fields1.clear();
                fields2.clear();
                        
                if ( new2old_illumina && !old_style_illumina_flag ) //if convert to old-style illumina headers is true and not old illumina files.
                {
                    split_str( line1, fields1, " " );
                    split_str( fields1[0], fields2, ":" );
                    line1 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) );
                            
                    fields1.clear();
                    fields2.clear();
                            
                    split_str( line2, fields1, " " );
                    split_str( fields1[0], fields2, ":" );
                            
                    line2 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) );//+ " (" + line2 + ")");
                            
                    fields1.clear();
                    fields2.clear();
                }
                        
                record_block1.push_back(line1); 
                record_block2.push_back(line2);
                        
                ii++;
                continue;
            }
            /*DNA string*/
            if(ii==1) 
            {
                record_block1.push_back(line1); /*DNA string*/
                record_block2.push_back(line2);
                ii++;
                continue;
            }
            /*a symbol "+"*/
            if(ii==2) 
            {
                record_block1.push_back(line1);
                record_block2.push_back(line2);
                ii++;
                continue;
            }
            if(ii==3) 
            {
                ii=0;
           
                Read *read1 = new Read();
                read1->illumina_readID = record_block1[0];
                read1->initial_length = record_block1[1].length();
                read1->read = record_block1[1];
                read1->illumina_quality_string = line1;
                        
                Read *read2 = new Read();
                read2->illumina_readID = record_block2[0];
                read2->initial_length = record_block2[1].length();
                read2->read = record_block2[1];
                read2->illumina_quality_string = line2;
                
                //Serial realization - useful for debugging if something does not work as expected
                for(int k =0; k<dict_holder.size(); k++) {
                        struct LibHitData match = CheckForLib2(&read1->read,&(dict_holder[k]));
                        if(match.start_pos != -1) 
                        {
                                string _str = match.lib_id + "\t" + int2str(match.start_pos) + "\t" + int2str(match.end_pos) + "\t" + read1->illumina_readID + "\t" + read1->read + "\t" + read1->illumina_quality_string + "\t" + read2->illumina_readID + "\t" + read2->read + "\t" + read2->illumina_quality_string + "\n";
                                fputs((char*)_str.c_str(), rep_file);
                  
                                WritePEFile(pe_output_file1, read1);
                                WritePEFile(pe_output_file2, read2);
                        }
                        else
                        {
                                match = CheckForLib2(&read2->read,&(dict_holder[k]));
                                if(match.start_pos != -1) 
                                {
                                        string _str = match.lib_id + "\t" + int2str(match.start_pos) + "\t" + int2str(match.end_pos) + "\t" + read1->illumina_readID + "\t" + read1->read + "\t" + read1->illumina_quality_string + "\t" + read2->illumina_readID + "\t" + read2->read + "\t" + read2->illumina_quality_string + "\n";
                                        fputs((char*)_str.c_str(), rep_file);
                        
                                        //Construct FASTQ entry
                                        WritePEFile(pe_output_file1, read1);
                                        WritePEFile(pe_output_file2, read2);
                                }
                        }
                }
                record_block1.clear();
                read1->illumina_readID.clear(); 
                read1->illumina_quality_string.clear();
                read1->read.clear();
          
                record_block2.clear();
                read2->illumina_readID.clear(); 
                read2->illumina_quality_string.clear();
                read2->read.clear();
          
                delete read1;
                delete read2;
                
            }
        }
        in1.close();
        in2.close();
    }
    
    
   
}


void IlluminaDynamicSE()
{
    fstream rep_file;
    fstream se_output_file;
    rep_file.open(rep_file_name.c_str(),ios::out);
    
    vector<string> record_block;
    
    for(int jj=0; jj<(int)se_names.size(); ++jj)
    {
        int ii = 0;
        
        std::string line;
        igzstream in(se_names[jj]);
        
        cout << "Processing file: " << se_names[jj] << "\n";
        
        while ( getline(in,line) )
        {
            /*Read ID*/
            if(ii==0) 
            {
                //Check for order
                vector <string> fields1, fields2;
                        
                if ( new2old_illumina && !old_style_illumina_flag ) //if convert to old-style illumina headers is true and not old illumina files.
                {
                    split_str( line, fields1, " " );
                    split_str( fields1[0], fields2, ":" );
                    line = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) );
                            
                    fields1.clear();
                    fields2.clear();
                }
                        
                record_block.push_back(line); 
                        
                ii++;
                continue;
            }
            /*DNA string*/
            if(ii==1) 
            {
                record_block.push_back(line); /*DNA string*/
                ii++;
                continue;
            }
            /*a symbol "+"*/
            if(ii==2) 
            {
                record_block.push_back(line);
                ii++;
                continue;
            }
            if(ii==3) 
            {
                ii=0;
           
                Read *read = new Read();
                read->illumina_readID = record_block[0];
                read->initial_length = record_block[1].length();
                read->read = record_block[1];
                read->illumina_quality_string = line;
                
                //Serial realization - useful for debugging if something does not work as expected
                for(int k =0; k<dict_holder.size(); k++) {
                        struct LibHitData match = CheckForLib2(&read->read,&dict_holder[k]);
                        if(match.start_pos != -1) 
                        {
                                rep_file << match.lib_id << "\t" << match.start_pos << "\t" << match.end_pos << "\t" << read->illumina_readID << "\t" << read->read << "\t" << read->illumina_quality_string << "\n";
                    
                                if(!se_output_file.is_open()) 
                                {
                                        se_output_file.open( se_output_filename.c_str(), ios::out );
                                }
                  
                                //Create FASTQ entry
                                WritePEFile(se_output_file, read);
                        } 
                }
                
                record_block.clear();
                read->illumina_readID.clear(); 
                read->illumina_quality_string.clear();
                read->read.clear();
                
                delete read;
                
            }
        }
        in.close();
        
    }
    
    
    cout << "====================Done====================\n";  
    
    
    rep_file.close();
}

