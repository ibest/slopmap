/* 
 * File:   util.h
 * Author: ilya
 *
 * Created on 31 Май 2012 г., 10:27
 */

#ifndef UTIL_H
#define	UTIL_H

#include <stdio.h>
#include <iostream>
#include "sparsehash-2.0.2/src/sparsehash/dense_hash_map"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <sys/stat.h>
#include <errno.h>
#include <sys/types.h>
#include <dirent.h>
#include <map>
#include <math.h>
 #include <stdint.h>

/* Numerical values for bases. */
#define T_BASE_VAL 0
#define U_BASE_VAL 0
#define C_BASE_VAL 1
#define A_BASE_VAL 2
#define G_BASE_VAL 3
#define N_BASE_VAL 4   /* Used in 1/2 byte representation. */

using namespace std;
//using namespace boost;
using google::dense_hash_map;
void stoupper(std::string& s);
char* itoa(int value, char* result, int base);
string MakeSeqComplement(string init_str);
string MakeRevComplement(string init_str);
//int GetRandomInt (int from, int to);
//double GetRandomDouble (double from, double to);
void split_str(const string& str, vector<string>& tokens, const string& delimiters);
void TrimNs(string &read);
void TrimNs2(string &read);
string GenNs(int num, char* letter);
int min4( int x, int y, int z, int k );
int max4( int x, int y, int z, int k );
int max3( int x, int y, int z );
int min3( int x, int y, int z );
bool exists(char* filePath);
double GetAvg( double past_avg, long n, int cur_diff );
int MakeDirectory(string path_to_create);
void GetDirectories(std::vector<string> &out, char *directory);
vector< vector<string> > GetPEfilenames(string prefix1, string prefix2, char *directory);
string i2str(int value, char* result, int base);
string double2str(double num);
string int2str(int num);
uint32_t dna_number(string s);
extern dense_hash_map<char, int > dnaDict;
#endif	/* UTIL_H */

