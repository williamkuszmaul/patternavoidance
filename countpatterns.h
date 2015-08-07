#ifndef _COUNTPATTERNS_H 
#define _COUNTPATTERNS_H // To avoid header being included twice in complilation process.

#include <assert.h>
#include <string.h>
//#include <iostream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <time.h> 
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <stdint.h>
#include <unordered_map>
#include <queue>
#include "hashdb.h"
#include "hashmap.h"
#include <sys/time.h>
#include "fastavoidance.h"
using namespace std;


void countpatterns(string patternlist, int maxpermsize, vector < vector <int> > & tally, vector < vector < int > > &completelist, bool verbose);

#endif 
