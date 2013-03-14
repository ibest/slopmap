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
#include <map>
#include <vector>
#include <pthread.h>
#include "util.h"
#include <streambuf>
#include <exception>
#include <pthread.h>
#include <math.h>
#include "Read.h"

using namespace std;

typedef struct  {
    //string seq_id; //seq id related to k_mer in the screening file
    string seq_id;
    int pos; //position of the k_mer in the screening file
} k_mer_struct;


/*Structure that hold seq_id and position of k_mer*/
typedef struct {
    vector <k_mer_struct> kmers;
    long pos;
    string k_mer_string;
    string string_to_align;
    
} HitData;


typedef struct {
    int pos;
    int flag;
    int i;
} DistanceStruct;

/*Extern variables*/
extern map<string, vector<k_mer_struct> > LibDict;
extern map<string, vector<k_mer_struct> >::iterator it_LibDict;
extern short KMER_SIZE;
extern short DISTANCE;
extern short NUM_CONSEQUITIVE_HITS;
extern short NUM_HITS;
extern bool mode_flag;


string CheckForLib(string seq);

#endif	/* KMERROUTINE_H */

