/* 
 * File:   KMerRoutine.h
 * Author: ilya
 *
 * Created on 7 Август 2012 г., 23:10
 */

#ifndef KMERROUTINE_H
#define	KMERROUTINE_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <pthread.h>
#include "util.h"
#include <streambuf>
#include <exception>
#include <pthread.h>
#include <signal.h>
#include <sys/select.h>
#include <math.h>
#include "Read.h"
#include <algorithm>    // std::sort
#include "dnautil.h"
#include "gzstream.h"
#include "sff.h"

/* D E F I N E S *************************************************************/
#define PRG_NAME "sff2fastq"
#define SFF_FILE_VERSION "0.8.8"


using namespace std;

typedef struct  {
    //string seq_id; //seq id related to k_mer in the screening file
    string seq_id;
    long pos; //position of the k_mer in the screening file
} k_mer_struct;


/*Structure that hold seq_id and position of k_mer*/
typedef struct {
    vector <k_mer_struct> kmers;
    long pos;
    string k_mer_string;
    string string_to_align;
    
} HitData;

struct LibHitData {
    long start_pos;
    long end_pos;
    string lib_id;
};


/*Extern variables*/
extern map<UBYTE, vector<k_mer_struct> > LibDict;
extern map<UBYTE, vector<k_mer_struct> >::iterator it_LibDict;
extern map<string, int > LibDictId; //Represent a pair: <LibId,length of the lib>
extern map<string, int >::iterator it_LibDictId;

extern short KMER_SIZE;
extern short DISTANCE;
//extern short NUM_CONSEQUITIVE_HITS;
//extern short NUM_HITS;
extern bool mode_flag;
extern double similarity_threshold;

extern char* illumina_file_name_R1;// = "";
extern char* illumina_file_name_R2;// = "";
extern char* illumina_file_name_se;
extern string pe_output_filename1;
extern string pe_output_filename2;
extern string se_output_filename;
extern string fastq_output_filename;
extern vector<char*> pe1_names, pe2_names, se_names, roche_names;
extern string output_prefix;
/*Report files*/
extern string rep_file_name;

extern FILE *rep_file;
extern fstream pe_output_file1, pe_output_file2;
    

extern bool old_style_illumina_flag;
extern int phred_coeff_illumina; //by default assume new illumina (1.8)
extern bool i64_flag;
extern bool new2old_illumina;

struct LibHitData CheckForLib(string seq);
struct LibHitData CheckForLib2(string *seq, map<UBYTE, vector<k_mer_struct> > *LibDict);

extern vector<map<UBYTE, vector<k_mer_struct> > > dict_holder;

void IlluminaDynamicSE();
void IlluminaDynamic();
void Roche454Dynamic();

void WriteFastqFile(fstream &fastq_output_file, string readID, char* bases, uint8_t* quality);

#endif	/* KMERROUTINE_H */

