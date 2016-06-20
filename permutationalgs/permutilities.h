#ifndef _UTILITIES_H 
#define _UTILITIES_H // To avoid header being included twice in complilation process.

#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>      /* printf, scanf, puts, NULL */
#include <cstdlib>
#include <cstring>
#include <ctime> 
#include <iostream>
#include <fstream>
#include <queue>
#include <sys/time.h>
#include <vector>
#include "perm.h"
#include "hashdb.h"
using namespace std;


// For patterns in S_{<10} can use like this:
//  hashdb patternset(1<<3); 
//  string permlist = "3124 4123 3142 4132";
//  makepatterns(permlist, patternset);
// Fills patternset with permutations in permlist, each as a 64 bit integer
void makepatterns(string permlist, hashdb &patternset, int &maxpatternsize);

// build prefix table containing all complements of normalizations of prefixes of every perm in permset
void addprefixes(const hashdb &permset, hashdb &table);

#endif 
