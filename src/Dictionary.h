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
//#include <sparsehash/sparse_hash_map>
//#include <sparsehash/sparse_hash_set>
//#include <sparsehash/dense_hash_set>
#include "Read.h"
#include "KMerRoutine.h"


using namespace std;

/*Extern variables*/
extern map<string, vector<k_mer_struct> > LibDict;
extern map<string, vector<k_mer_struct> >::iterator it_LibDict;
extern map<string, int > LibDictId; //Represent a pair: <LibId,length of the lib>
extern short KMER_SIZE;
 
/*Builds a new Contaminants Dictionary. Here it assumes that frequency is 1*/
int BuildLibDictionary(string filename);

void PutLibKmer(string str, string contaminant_id);

#endif	/* DICTIONARY_H */

