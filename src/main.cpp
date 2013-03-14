#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <streambuf>
#include <exception>
#include <pthread.h>
#include <bits/basic_string.h>
#include "timer.h"
#include "util.h"
#include "Read.h"
#include "Dictionary.h"
#include "KMerRoutine.h"
#include "gzstream.h"
#include "sff.h"
#include "ascii.h"
#include <stdlib.h>


using namespace std;

/* D E F I N E S *************************************************************/
#define VERSION "0.8.0"
#define PRG_NAME "sff2fastq"


string version = "0.3.3 (2013-03-13)"; 

/*Computational parameters (default)*/
short KMER_SIZE = 15;
short DISTANCE = 5;
short NUM_CONSEQUITIVE_HITS = 3;
short NUM_HITS = 10;
bool mode_flag = false; //If this flag is true, the program will check NUM_HITS instead of NUM_CONSEQUITIVE_HITS.

map<string, vector<k_mer_struct> > LibDict;
map<string, vector<k_mer_struct> >::iterator it_LibDict;

char *lib_filename;
bool lib_flag = false;

/*Illumina*/
bool illumina_flag = false;
bool illumina_se_flag = false;
char* illumina_file_name_R1;// = "";
char* illumina_file_name_R2;// = "";
char* illumina_file_name_se;
vector<char*> pe1_names, pe2_names, se_names, roche_names;

/*Roche 454*/
bool roche_flag = false; 


bool old_style_illumina_flag = false;
int phred_coeff_illumina = 33; //by default assume new illumina (1.8)
bool i64_flag = false;
bool new2old_illumina = false;

string output_prefix;

/*Report files*/
string rep_file_name;

void Roche454Dynamic();
void IlluminaDynamic();
void IlluminaDynamicSE();


int main(int argc, char *argv[]) 
{
    double start, finish, elapsed;
    GET_TIME(start);
    
    cout << "Version: " << version << endl;
    /*******************************************/
    /* Parse command line arguments */
    /*******************************************/
    if(argv[1] == NULL) {
        //PrintHelp();
        return 0;
    }
    
    if( (string(argv[1]) == "-help") || (string(argv[1]) == "--help") || (string(argv[1]) == "-?") ) {
        //PrintHelp();
        return 0;
    }
        
    for (int i=1; i<argc; i++) 
    {
        if( string(argv[i]) == "--version" ) 
        {
           cout << "Version: " << version << endl;
           exit(1);
        }
        
        if( string(argv[i]) == "-k" ) 
        {
           KMER_SIZE = atoi(argv[++i]);
           continue;
        }
        if( string(argv[i]) == "-d" ) 
        {
           DISTANCE = atoi(argv[++i]);
           continue;
        }
        if( string(argv[i]) == "--mode2" ) 
        {
           mode_flag = true;
           continue;
        }
        if( string(argv[i]) == "-nch" ) 
        {
           NUM_CONSEQUITIVE_HITS = atoi(argv[++i]);
           continue;
        }
        if( string(argv[i]) == "-nh" ) 
        {
           NUM_HITS = atoi(argv[++i]);
           continue;
        }
        if(string(argv[i]) == "-l" )
        { 
           if ((i+1)<argc ) 
           {
             lib_filename = argv[++i]; /*File with contaminants given*/
             lib_flag = true;
           }
           continue;
        }
        if( string(argv[i]) == "-1" )
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_flag = true;
              illumina_file_name_R1 = argv[++i];
              pe1_names.push_back(illumina_file_name_R1);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  pe1_names.push_back(argv[i+jj+1]);
                  //printf("%s\n", argv[i+jj+1]);
                  jj+=1;
              }
              
           }
           
           continue;
        }
        if( string(argv[i]) == "-2" )
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_flag = true;
              illumina_file_name_R2 = argv[++i];
              pe2_names.push_back(illumina_file_name_R2);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  pe2_names.push_back(argv[i+1+jj]);
                  //printf("%s\n", argv[i+1+jj]);
                  jj+=1;
              }
              
           }
           
           continue;
        } 
        if( string(argv[i]) == "-U" ) //single-end file mode
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_se_flag = true;
              illumina_flag = true;
              illumina_file_name_se = argv[++i];
              se_names.push_back(illumina_file_name_se);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  se_names.push_back(argv[i+1+jj]);
                  //printf("%s\n", argv[i+1+jj]);
                  jj+=1;
              }
           }
           
           continue;
        }
        if( string(argv[i]) == "-454" ) 
        {
            if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
            {
              roche_flag = true; 
              
              roche_names.push_back(argv[++i]);
              
              int jj=0;
              while( ( (i+1+jj)<argc ) && (argv[i+1+jj][0] != '-') )
              {
                  roche_names.push_back(argv[i+1+jj]);
                  jj+=1;
              }
              
            }
            
            continue;
        }
        if(string(argv[i]) == "-o" ) 
        {
           if(argv[i+1] == NULL) 
           {
              cout << "Error: you have to specify the output prefix of output files." << endl;
              return -1;
           } else
           {
               output_prefix = string(argv[++i]);
           }
        }
    }
    
    if(output_prefix == "")
    {
        cout << "No output prefix found.\n";
        return -1;
    }
    
    //Check if input files exist
    if (illumina_flag)
    {
        if(illumina_se_flag)
        {
            for(int i=0; i<(int)se_names.size(); ++i)
            {
                if ( !exists( se_names[i] ) )
                {
                  cout<< "Error: file " <<  se_names[i] << " does not exist\n";
                  return 0;
                }
                
                //Test is the files provided are old-style illumina
                std::string line;
                igzstream in(se_names[i]); 
                getline(in,line);
                vector <string> fields;
                split_str( line, fields, ":" );
                if (fields.size() == 5)
                {
                    //Old headers == True
                    old_style_illumina_flag = true;
                    
                } else //if (fields.size() == 7)
                {
                    //Old headers == False
                    old_style_illumina_flag = false;
                }
                
            }
        } else 
        {
                if(pe1_names.size() != pe2_names.size())
                {
                        cout<< "Error: numbers of PE1 files and PE2 files do not match!\n";
                        return 0;
                } else
                {
                        for(int i=0; i<(int)pe1_names.size(); ++i)
                        {
        
                                if ( !exists( pe1_names[i] ) )
                                {
                                        cout<< "Error: file " <<  pe1_names[i] << " does not exist\n";
                                        return 0;
                                }
                                if (!exists( pe2_names[i] ) )
                                {
                                        cout<< "Error: file " <<  pe2_names[i] << " does not exist\n";
                                        return 0;
                                }
                                
                                
                                //Test is the files provided are old-style illumina
                                std::string line1, line2;
                                igzstream in1(pe1_names[i]); igzstream in2(pe2_names[i]); 
                                getline(in1,line1); getline(in2,line2);
                                vector <string> fields1, fields2;
                                split_str( line1, fields1, ":" ); split_str( line2, fields2, ":" );
                                if ( (fields1.size() == 5) && (fields2.size() == 5) )
                                {
                                        //Old headers == True
                                        old_style_illumina_flag = true;
                                        
                                } else if ( (fields1.size() == 10) && (fields2.size() == 10) )
                                {
                                        //Old headers == False
                                        old_style_illumina_flag = false;
                                } else if ( fields1.size() != fields2.size() )
                                {
                                    cout << "Error: impossible to have both old & new illumina as paired-end reads" << endl;
                                    return -1;
                                } else 
                                {
                                    cout << "Warning: unknown header format in file: " << pe1_names[i] << ", " << pe2_names[i] << endl;
                                     old_style_illumina_flag = false;
                                }
                        }
                        
                }
        }
        
        
        if (i64_flag == true)
           phred_coeff_illumina = 64;
        
    } 
    
    if(roche_flag)
    {
        for(int i=0; i<(int)roche_names.size(); ++i)
        {
            if (!exists( roche_names[i] ) )
            {
               cout<< "Error: file " <<  roche_names[i] << " does not exist\n";
               return 0;
            }
        }
    }
    
    if(!illumina_flag && !roche_flag)
    {
        printf("Error! You have to specify R1 or R2 files or both or Roche 454 file (sff or fastq)\n");
        return(-1);
    }
    
    if(lib_flag)
    {
        if (!exists( lib_filename ) )
        {
            cout<< "Error: library file provided " <<  lib_filename << " does not exist\n";
            return -1;
        }
    } 
    else
    {
        printf("Error! You have to specify the library file.\n");
        return(-1);
    }
    
    /*Building dictionary*/
    BuildLibDictionary(lib_filename);
    
    rep_file_name = output_prefix + ".txt" ;
    
    if(!illumina_se_flag)
    {
        IlluminaDynamic();
    }
    else
    {
        IlluminaDynamicSE();
    }
    
    if(roche_flag)
    {
        Roche454Dynamic();
    }
    
    GET_TIME(finish);
    elapsed = finish - start;
    
    printf("Elapsed time = %e seconds\n", elapsed);
    printf("Program finished.\n");
}

void Roche454Dynamic()
{
    fstream rep_file;
    rep_file.open(rep_file_name.c_str(),ios::out);
    
    for(int i=0; i<(int)roche_names.size(); ++i)
    {
            cout << "Parsing file: " << roche_names[i] << "..." << endl;
        
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
                verify_sff_common_header((char*)PRG_NAME, (char*)VERSION, &h);

                int left_clip = 0, right_clip = 0, nbases = 0;
                char *bases;
                uint8_t *qualily;
                
                int numreads = (int) h.nreads;
                
                for (int i = 0; i < numreads; i++) 
                { 
                    read_sff_read_header(sff_fp, &rh);
                    read_sff_read_data(sff_fp, &rd, h.flow_len, rh.nbases);

                    get_clip_values(rh, 0, &left_clip, &right_clip);
                    nbases = right_clip - left_clip;

                    // create bases string 
                    bases = get_read_bases(rd, left_clip, right_clip);

                    // create quality array 
                    qualily = get_read_quality_values(rd, left_clip, right_clip);

                    for (int j = 0; j < nbases; j++) 
                    {
                        qualily[j] = (qualily[j] <= 93 ? qualily[j] : 93) + 33;
                    }
                    
                    string lib_id = CheckForLib(string(bases));
                    if(lib_id != "") 
                    {
                        rep_file << lib_id << "\t" << rh.name << "\t" << bases << "\t" << qualily << "\n";
                    }
                    
                    
                    free(bases);
                    free(qualily);
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
                       quality = line;
                       
                       string lib_id = CheckForLib(string(bases));
                       if(lib_id.length() > 0) 
                       {
                          rep_file << lib_id << "\t" << readID << "\t" << bases << "\t" << quality << "\n";
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
    FILE *rep_file;

    if ( (rep_file = fopen((char*)rep_file_name.c_str(), "w")) == NULL ) 
    {
        fprintf(stderr,"[err] Could not open file '%s' for reading.\n", (char*)rep_file_name.c_str());
        exit(1);
    }
    
    vector<string> record_block1, record_block2;
    
    for(int jj=0; jj<(int)pe1_names.size(); ++jj)
    {
        int ii = 0;
        
        std::string line1, line2;
        igzstream in1(/*fastq_file1*/pe1_names[jj]); //for R1
        igzstream in2(/*fastq_file2*/pe2_names[jj]); //for R2
        
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
                    line1 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) + " (" + line1 + ")");
                            
                    fields1.clear();
                    fields2.clear();
                            
                    split_str( line2, fields1, " " );
                    split_str( fields1[0], fields2, ":" );
                            
                    line2 = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) + " (" + line2 + ")");
                            
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
                string lib_id = CheckForLib(read1->read);
                if(lib_id.length() > 0) 
                {
                    string _str = lib_id + "\t" + read1->illumina_readID + "\t" + read1->read + "\t" + read1->illumina_quality_string + "\t" + read2->illumina_readID + "\t" + read2->read + "\t" + read2->illumina_quality_string + "\n";
                    fputs((char*)_str.c_str(), rep_file);
                    //rep_file << lib_id << "\t" << read1->illumina_readID << "\t" << read1->read << "\t" << read1->illumina_quality_string << "\t" << read2->illumina_readID << "\t" << read2->read << "\t" << read2->illumina_quality_string << "\n";
                } 
                else 
                {
                    lib_id = CheckForLib(read2->read);
                    if(lib_id.length() > 0) {
                        //rep_file << lib_id << "\t" << read1->illumina_readID << "\t" << read1->read << "\t" << read1->illumina_quality_string << "\t" << read2->illumina_readID << "\t" << read2->read << "\t" << read2->illumina_quality_string << "\n";
                        string _str = lib_id + "\t" + read1->illumina_readID + "\t" + read1->read + "\t" + read1->illumina_quality_string + "\t" + read2->illumina_readID + "\t" + read2->read + "\t" + read2->illumina_quality_string + "\n";
                        fputs((char*)_str.c_str(), rep_file);
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
    
    
    cout << "====================Done====================\n";  
    fclose(rep_file);
   
}


void IlluminaDynamicSE()
{
    fstream rep_file;
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
                    line = string(fields2[0] + "_" + fields2[2] + ":" + fields2[3] + ":" + fields2[4] + ":" + fields2[5] + ":" + fields2[6] + "#0/" + fields1[1].substr(0,1) + " (" + line + ")");
                            
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
                string lib_id = CheckForLib(read->read);
                if(lib_id.length() > 0) 
                {
                    rep_file << lib_id << "\t" << read->illumina_readID << "\t" << read->read << "\t" << read->illumina_quality_string << "\n";
                    if (rep_file.fail())
                        {
                                std::cout << "failed to append to file!\n";
                        
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
