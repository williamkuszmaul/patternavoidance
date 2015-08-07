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



typedef unsigned long long timestamp_t;

// copied from http://stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c
static timestamp_t
    get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}


// indexing starts at 0
inline uint64_t getdigit(uint64_t perm, int index) {
  return (perm >> (index * 4))  & 15L; // grab digit
}

inline uint64_t setdigit(uint64_t perm, int index, uint64_t newdigit) {
  return (perm & ~(15L << (index * 4))) |  (newdigit << (index * 4)); // clear digit and then rewrite its value
}

void displayperm(uint64_t perm);

void displayperm(uint64_t perm, int size);

// kills digit in position index // DOES NOT NORMALIZE PERMUTATION AFTERWARDS
inline uint64_t killpos(uint64_t perm, int index) {
  uint64_t bottom = perm & ((1L << (4 * index)) - 1);
  uint64_t top = perm & ((~ 0L) - ((1L << (4 * index + 4)) - 1));
  if (index == 15) return bottom; // top is ill-defined in this case
  return bottom + (top >> 4); 
}

inline uint64_t addpos(uint64_t perm, int index) {
  uint64_t bottom = perm & ((1L << (4 * index)) - 1);
  uint64_t top = perm & ((~ 0L) - ((1L << (4 * index)) - 1));
  return bottom + (top << 4);
}


inline uint64_t getinverse(uint64_t perm, int length) {
  uint64_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, getdigit(perm, i), i);
  }
  return answer;
}


uint64_t makepatterns(string permlist, hashdb &patternset, int &maxpatternsize);

void buildavoiders(const hashdb &patternset, int maxavoidsize, int maxsize, vector < vector < uint64_t > > &avoidervectorm, uint64_t plannedavoidsetsize);

void countavoiders(const hashdb &patternset, int maxavoidsize, int maxsize, vector < uint64_t > &numavoiders, uint64_t plannedavoidsetsize);

void buildavoidersfrompatternlist(string patternlist, int maxpermsize, vector < vector < uint64_t > > &avoidervector);

void countavoidersfrompatternlist(string patternlist, int maxpermsize, vector < uint64_t > &numavoiders);
//void countavoidersv2(const hashdb &patternset, int maxavoidsize, int maxsize, vector < uint64_t > &numavoiders, uint64_t plannedavoidsetsize);

void countavoidersfromfile(ifstream &infile, ofstream &outfile, int maxpermsize, bool verbose);

#endif 
