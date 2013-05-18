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
#include "sparsehash-2.0.2/src/google/dense_hash_map"
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

using namespace std;
using google::dense_hash_map;      // namespace where class lives by default
using tr1::hash;  // or __gnu_cxx::hash, or maybe tr1::hash, depending on your OS

/*Extern variables*/
extern dense_hash_map<string, vector<k_mer_struct> > LibDict;
extern dense_hash_map<string, vector<k_mer_struct> >::iterator it_LibDict;
extern dense_hash_map<string, int > LibDictId; //Represent a pair: <LibId,length of the lib>
extern short KMER_SIZE;
 
/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildLibDictionary(char* filename);

void PutLibKmer(string str, string contaminant_id);

#endif	/* DICTIONARY_H */

