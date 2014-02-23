/* 
 * File:   Dictionary.h
 * Author: ilya
 *
 * Created on 7 Август 2012 г., 0:48
 */

#ifndef DICTIONARY_H
#define	DICTIONARY_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <pthread.h>
#include "util.h"
#include <streambuf>
#include <exception>
#include <pthread.h>
#include <math.h>
#include "Read.h"
#include "KMerRoutine.h"
#include "dnautil.h"
 #include <stdint.h>


using namespace std;

/*Extern variables*/
extern map<UBYTE, vector<k_mer_struct> > LibDict;
extern map<UBYTE, vector<k_mer_struct> >::iterator it_LibDict;
extern map<string, int > LibDictId; //Represent a pair: <LibId,length of the lib>
extern short KMER_SIZE;
extern bool illumina_pe_flag;
extern bool illumina_se_flag;
extern bool roche_flag;

extern FILE *rep_file;
extern fstream pe_output_file1, pe_output_file2;

extern vector<map<UBYTE, vector<k_mer_struct> > > dict_holder;


/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildLibDictionary(char* filename);
int BuildLibDictionary2(char* filename);
struct twoBit *twoBitFromDnaSeq(struct dnaSeq *seq, string lib_id);
void PutLibKmer(string str, string contaminant_id);

#endif	/* DICTIONARY_H */

