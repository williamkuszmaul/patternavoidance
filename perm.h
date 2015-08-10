// A permutation is represented as a 64 bit integer, with the ith digit in the ith half-byte.
// Only permutations of size <= 16 can be represented this way
// This file contains basic functions for this permuation representation
// It also, kind of randomly, contains a function for getting the precise time of day

#ifndef _PERM_H 
#define _PERM_H // To avoid header being included twice in complilation process.

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
#include "hashdb.h"
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

// indexing starts at 0 for both position and value in permutation
inline uint64_t getdigit(uint64_t perm, int index) {
  return (perm >> (index * 4))  & 15L; // grab digit
}

// indexing starts at 0 for both position and value in permutation
inline uint64_t setdigit(uint64_t perm, int index, uint64_t newdigit) {
  return (perm & ~(15L << (index * 4))) |  (newdigit << (index * 4)); // clear digit and then rewrite its value
}

inline void displayperm(uint64_t perm) {
  for (int i = 0; i < 16; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

inline void displayperm(uint64_t perm, int size) {
  for (int i = 0; i < size; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}


// kills digit in position index (with indexing startign at zero)
// DOES NOT NORMALIZE PERMUTATION AFTERWARDS
inline uint64_t killpos(uint64_t perm, int index) {
  uint64_t bottom = perm & ((1L << (4 * index)) - 1);
  uint64_t top = perm & ((~ 0L) - ((1L << (4 * index + 4)) - 1));
  if (index == 15) return bottom; // top is ill-defined in this case
  return bottom + (top >> 4); 
}

// inserts a blank in position index
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


inline uint64_t getreverse(uint64_t perm, int length) {
  uint64_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, i, getdigit(perm, length - i - 1));
  }
  return answer;
}

inline uint64_t getcomplement(uint64_t perm, int length) {
  uint64_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, i, length - getdigit(perm, i) - 1);
  }
  return answer;
}

// For patterns in S_{<10} can use like this:
//  hashdb patternset(1<<3); 
//  string permlist = "3124 4123 3142 4132";
//  makepatterns(permlist, patternset);
// Fills patternset with permutations in permlist, each as a 64 bit integer
void makepatterns(string permlist, hashdb &patternset, int &maxpatternsize);

// gives permuation corresponding with a string
inline uint64_t stringtoperm(string perm) {
  uint64_t answer = 0;
  for (int i = 0; i < perm.size(); i++) answer = setdigit(answer, i, (unsigned long long)(perm[i] - '0'));
  return answer;
}

// A fast map bijection from S_length to \mathbb{Z}_{length!}
unsigned long long permtonum(uint64_t perm, int length);




#endif 
