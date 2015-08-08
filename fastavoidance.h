// fastavoidance.h


#ifndef _FASTAVOIDANCE_H 
#define _FASTAVOIDANCE_H // To avoid header being included twice in complilation process.


#include <assert.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <time.h> 
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <stdint.h>
// #include <unordered_set>
#include <queue>
#include "hashdb.h"
#include <sys/time.h>
using namespace std;



void buildavoiders(const hashdb &patternset, int maxavoidsize, int maxsize,  vector < vector < uint64_t > > &avoidervector, vector < uint64_t > &numavoiders, bool justcount, uint64_t plannedavoidsetsize);

void buildavoidersfrompatternlist(string patternlist, int maxpermsize, vector < vector < uint64_t > > &avoidervector);

void countavoidersfrompatternlist(string patternlist, int maxpermsize, vector < uint64_t > &numavoiders);
//void countavoidersv2(const hashdb &patternset, int maxavoidsize, int maxsize, vector < uint64_t > &numavoiders, uint64_t plannedavoidsetsize);

void countavoidersfromfile(ifstream &infile, ofstream &outfile, int maxpermsize, bool verbose);

#endif 
