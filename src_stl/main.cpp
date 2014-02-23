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
#include "ascii.h"
#include <stdlib.h>
#include <stdint.h>

using namespace std;

string version = "2.0.2 (2014-02-23)"; 

/*Computational parameters (default)*/
short KMER_SIZE = 15;
short DISTANCE = 1;
bool mode_flag = false; //If this flag is true, the program will check NUM_HITS instead of NUM_CONSEQUITIVE_HITS.

char *lib_filename;
bool lib_flag = false;

/*Illumina*/
bool illumina_pe_flag = false;
bool illumina_se_flag = false;
char* illumina_file_name_R1;// = "";
char* illumina_file_name_R2;// = "";
char* illumina_file_name_se;
string pe_output_filename1;
string pe_output_filename2;
string se_output_filename;
string fastq_output_filename;

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

double similarity_threshold = 0.75; //75% of similarity

void Roche454Dynamic();
void IlluminaDynamic();
void IlluminaDynamicSE();
void WritePEFile(fstream &pe_output_file, Read *read);
void Roche454Dynamic2();

FILE *rep_file;
fstream pe_output_file1, pe_output_file2;
    
//vector<dense_hash_map<UBYTE, vector<k_mer_struct> > > dict_holder;
vector<map<UBYTE, vector<k_mer_struct> > > dict_holder;

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
        if( string(argv[i]) == "-t" ) 
        {
           similarity_threshold = atof(argv[++i]);
           if(similarity_threshold > 1) 
           {
               cout << "Threshold must be <= 1\n";
               return -1;
           }
           continue;
        }
        if( string(argv[i]) == "-k" ) 
        {
           KMER_SIZE = atoi(argv[++i]); //kmer length for sampling both the dictionary and the read
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
        if(string(argv[i]) == "-l" )
        { 
           if ((i+1)<argc ) 
           {
             lib_filename = argv[++i]; /*File with contaminants given*/
             lib_flag = true;
           }
           continue;
        }
        if( string(argv[i]) == "-1" ) //Read #1, Illumina
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_pe_flag = true;
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
        if( string(argv[i]) == "-2" ) //Read #2, Illumina
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_pe_flag = true;
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
        if( string(argv[i]) == "-U" ) //single-end read mode (Illumina)
        {
           if ( ( (i+1)<argc ) && (argv[i+1][0] != '-') ) 
           {
              illumina_se_flag = true;
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
        if( string(argv[i]) == "-454" ) //Roche 454 mode
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
    
    //************************VERIFICATION OF THE USER'S INPUTS*****************************
    
    if(output_prefix == "")
    {
        cout << "No output prefix found.\n";
        return -1;
    }
    
    
    if (i64_flag == true)
       phred_coeff_illumina = 64;
    
    //Check if input files exist
    if (illumina_pe_flag)
    {
        if(pe1_names.size() != pe2_names.size())
        {
            cout<< "Error: numbers of PE1 files and PE2 files do not match!\n";
            return 0;
        }
        
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
                        
       pe_output_filename1 =  output_prefix + "_PE1.fastq";
       pe_output_filename2 =  output_prefix + "_PE2.fastq";
          
    } else if(illumina_se_flag) {
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
            se_output_filename =  output_prefix + "_SE.fastq";
            
    } else if(roche_flag) {
        for(int i=0; i<(int)roche_names.size(); ++i)
        {
            if (!exists( roche_names[i] ) )
            {
               cout<< "Error: file " <<  roche_names[i] << " does not exist\n";
               return 0;
            }
        }
    }
    
    if(!illumina_pe_flag && !roche_flag)
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
    
    //Check if output path exist
    vector<string> t;
    split_str(output_prefix, t, "/");
    string t_prefix = "";
    if(output_prefix[0] == '/') t_prefix += '/';
    for(int ii=0; ii < (int)t.size()-1; ii++)
    {
        t_prefix += t[ii] + "/";
    }
    if ( (t_prefix != "") && (!exists( (char*)t_prefix.c_str() ) ) )
    {
         cout<< "Warning: path " <<  t_prefix << " in output prefix does not exist" << endl;
         cout << "Trying to create...\n";
         if ( MakeDirectory(t_prefix) == 0 )
         {
             cout << "Sucess!\n";
         } 
         else
         {
             cout << "Could not created following path: " << t_prefix << "\n";
             return 0;
         }
    }
    t.clear();
    t_prefix.clear();        
    
    
    
    rep_file_name = output_prefix + ".txt" ;
    fastq_output_filename = output_prefix + ".fastq" ;
    
    /*Building dictionary*/
    BuildLibDictionary2(lib_filename);
    
    
    GET_TIME(finish);
    elapsed = finish - start;
    
    printf("Elapsed time = %e seconds\n", elapsed);
    printf("Program finished.\n");
}

