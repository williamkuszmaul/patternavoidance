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
#include <boost/multiprecision/cpp_int.hpp>
#include "math.h"
#include <sys/time.h>
#include <vector>
using namespace std;

// Choices are 16, 25, or 42. All three choices work for avoidance
// code.  However, third choice is much slower than the others
#define MAXPERMSIZE 16

// permutations up to size 16
#if (MAXPERMSIZE == 16)
#define LETTERSIZE 4
typedef uint64_t perm_t;
#define numbits 64
#endif
// permutations up to size 25
#if (MAXPERMSIZE == 25)
#define LETTERSIZE 5
typedef boost::multiprecision::uint128_t perm_t;
#define numbits 128
#endif
// permutations up to size 42
#if (MAXPERMSIZE == 42)
#define LETTERSIZE 6
typedef boost::multiprecision::uint256_t perm_t;
#define numbits 256
#endif

#define LETTERFACE ((((perm_t) 1) << LETTERSIZE) - (perm_t)1)


//  NOTE: WE USE SIZEOF(PERM_T) IN WAYS THAT WILL NOT BE VALID IF THERE IS ANY EXTRA BAGGAGE ON PERM_T, BUT ONLY IN THIS FILE
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
inline uint64_t getdigit(perm_t perm, int index) {
  return (uint64_t)((perm >> (index * LETTERSIZE))  & LETTERFACE); // grab digit
}

// indexing starts at 0 for both position and value in permutation
inline perm_t setdigit(perm_t perm, int index, uint64_t newdigit) {
  return (perm & ~((perm_t)LETTERFACE << (index * LETTERSIZE))) |  ((perm_t)newdigit << (index * LETTERSIZE)); // clear digit and then rewrite its value
}

inline perm_t incrementdigit(perm_t perm, int index) {
  return (perm + ((perm_t)1 << (LETTERSIZE * index)));
}

inline perm_t decrementdigit(perm_t perm, int index) {
  return (perm - ((perm_t)1 << (LETTERSIZE * index)));
}



inline void displayperm(perm_t perm) {
  for (int i = 0; i < numbits / LETTERSIZE; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

inline void displayperm(perm_t perm, int size) {
  for (int i = 0; i < size; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}


// kills digit in position index (with indexing startign at zero)
// DOES NOT NORMALIZE PERMUTATION AFTERWARDS
inline perm_t killpos(perm_t perm, int index) {
  perm_t bottom = perm & (((perm_t)1 << (LETTERSIZE * index)) - 1);
  perm_t top = perm & ((~ (perm_t)0) - (((perm_t)1 << (LETTERSIZE * index + LETTERSIZE)) - 1));
  if ((LETTERSIZE * index + LETTERSIZE) == numbits) return bottom; // top is ill-defined in this case
  return bottom + (top >> LETTERSIZE); 
}

// inserts a blank in position index
inline perm_t addpos(perm_t perm, int index) {
  perm_t bottom = perm & (((perm_t)1 << (LETTERSIZE * index)) - 1);
  perm_t top = perm & ((~ (perm_t)0) - (((perm_t)1 << (LETTERSIZE * index)) - 1));
  return bottom + (top << LETTERSIZE);
}


inline perm_t getinverse(perm_t perm, int length) {
  perm_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, getdigit(perm, i), i);
  }
  return answer;
}


inline perm_t getreverse(perm_t perm, int length) {
  perm_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, i, getdigit(perm, length - i - 1));
  }
  return answer;
}

inline perm_t getcomplement(perm_t perm, int length) {
  perm_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, i, length - getdigit(perm, i) - 1);
  }
  return answer;
}

// gives permuation corresponding with a string
inline perm_t stringtoperm(string perm) {
  perm_t answer = 0;
  for (int i = 0; i < perm.size(); i++) answer = setdigit(answer, i, (unsigned long long)(perm[i] - '1'));
  return answer;
}

// A fast map bijection from S_length to \mathbb{Z}_{length!}
unsigned long long permtonum(perm_t perm, int length);

inline uint64_t getmaxdigit(perm_t perm) {
  uint64_t answer = 0;
  for (int i = 0; i < numbits / LETTERSIZE; i++) {
    answer = max(answer, getdigit(perm, i));
  }
  return answer;
}


// Input: perm, perm's inverse (which needs to be correct in position length - index), perm's length, index, answer = complement of normalization of (index)-prefix of perm, a bitmap named seenpos which should start off at zero for index = 0. 
// Output: bitmap is updated to keep track of the positions in perm of each letter from n - i to n. answer is updated to be complement of normalization of (index + 1)-prefix of perm
void extendnormalizetop(perm_t perm, perm_t inverse, int length, int index, perm_t &answer, uint32_t & seenpos);

// works if key is 64 bit integer
// requires maxsize is a power of two. Returns a number from zero to maxsize - 1.
inline unsigned long long hash_perm (perm_t key_in, uint64_t maxsize)
{
#if (numbits <= 64)
  uint64_t key = key_in;
#endif
#if (numbits == 128)
  perm_t bits = ((perm_t) 1 << 64) - 1;
  uint64_t key = (uint64_t) (bits & (key_in + (key_in >> 61))); // for 128 bit case
#endif
#if (numbits == 256)
  perm_t bits = ((perm_t) 1 << 64) - 1;
  uint64_t key = (uint64_t) (bits & (key_in + (key_in >> 61) + (key_in >> 126) + (key_in >> 189)));
#endif
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return ((uint64_t)key) & (maxsize-1); //assume maxsize is power of two
  //This hash function comes from http://www.concentric.net/~ttwang/tech/inthash.htm
}

inline perm_t make_key_nonzero (perm_t key) {
  return key + 1; // for all permutations we consider, adding 1 is guaranteed to result in a non-zero integer
}

inline perm_t revert_stored_key (perm_t key) {
  return key - 1;
}

inline bool not_zero_perm(perm_t key) {
#if (numbits == 256) // In this case, because of how the boost library works, being zero as an integer is not the same as having all zero bytes.
  for (int i = 0; i < sizeof(key) / 8; i++) {
    if (*(((uint64_t *)(&key)) + i) != 0L) {
	return true;
    }
  }
  return false;
#endif
  return key != (perm_t)0;
}

inline void assert_key_nonzero (perm_t key) {
  assert(not_zero_perm(key));
}

#endif 
