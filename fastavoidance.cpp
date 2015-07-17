// simplest make is g++ -O2 -std=c++11  fastavoidance.cpp  -g -o fastavoidance
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

#define ULL (uint64_t)
#define one ((uint64_t) 1)
#define zero ((uint64_t) 0)
#define bits4 ((uint64_t) 15)


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
  return (perm >> (index * 4))  & bits4; // grab digit
}

inline uint64_t setdigit(uint64_t perm, int index, uint64_t newdigit) {
  return (perm & ~(bits4 << (index * 4))) |  (newdigit << (index * 4)); // clear digit and then rewrite its value
}


void displayperm(uint64_t perm) {
  for (int i = 0; i < 16; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

void displayperm(uint64_t perm, int size) {
  for (int i = 0; i < size; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

// kills digit in position index // DOES NOT NORMALIZE PERMUTATION AFTERWARDS
inline uint64_t killpos(uint64_t perm, int index) {
  uint64_t bottom = perm & ((one << (4 * index)) - 1);
  uint64_t top = perm & ((~ zero) - ((one << (4 * index + 4)) - 1));
  return bottom + (top >> 4); 
}

inline uint64_t addpos(uint64_t perm, int index) {
  uint64_t bottom = perm & ((one << (4 * index)) - 1);
  uint64_t top = perm & ((~ zero) - ((one << (4 * index)) - 1));
  return bottom + (top << 4);
}


inline uint64_t getinverse(uint64_t perm, int length) {
  uint64_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, getdigit(perm, i), i);
  }
  return answer;
}

bool isavoider(uint64_t perm, int maxavoidsize, int length, hashdb &avoidset, hashdb &patternset) { 
  uint64_t inverse = getinverse(perm, length);
  if (patternset.contains(perm)) { // if is in set of bad patterns
      return false;
  }
  uint64_t currentperm = perm;
  if (length > 1) { // don't deal with permutations of size zero
    for (int i = length - 1; i >= 0 && i >= length - maxavoidsize - 1; i--) {
      if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
        currentperm = addpos(currentperm, getdigit(inverse, i + 1));
	currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
      }
      currentperm = killpos(currentperm, getdigit(inverse, i));
      //displayperm(currentperm);
      if (!avoidset.contains(currentperm)) {
	return false; // found a subword not avoiding the patterns
      }
    }
  }
  return true;
}

int countavoiders(hashdb &patternset, int maxavoidsize, int maxsize) {
  int numavoiders = 0;
  
  hashdb avoidset = hashdb(1<<20);
  uint64_t startperm = 0;
  avoidset.add(0); // identity in S_1
  
  std::queue<long long> avoiderstoextend;
  avoiderstoextend.push(startperm);

  int currentlength = 1; // maintain as length of next permutation to be popped from avoiderstoextend
  int numleftcurrentlength = 1; // number of permutations left in avoiderstoextend until we have to increment currentlength
  int numnextlength = 0; // number of permutations of size currentlength + 1 in avoiderstoextend

  while (avoiderstoextend.size() > 0) {
    if (numleftcurrentlength == 0) {
      numleftcurrentlength = numnextlength;
      numnextlength = 0;
      currentlength++;
    }
    uint64_t perm = avoiderstoextend.front();
    avoiderstoextend.pop();
    numleftcurrentlength--;
    if (currentlength >= maxsize) {
      break;
    }    
    for (int i = 0; i < currentlength + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentlength);
      if (isavoider(extendedperm, maxavoidsize, currentlength + 1, avoidset, patternset)) {
      	if (currentlength + 1 == maxsize) numavoiders++;
      	avoiderstoextend.push(extendedperm);
	numnextlength++;
	avoidset.add(extendedperm);
      }
    }
  }
  return numavoiders;
}




int main() {
  int maxpatternsize = 3;
  int permsize = 10;
  assert(permsize <= 16);
  uint64_t perm = 0;
  perm = setdigit(perm, 0, 0);
  perm = setdigit(perm, 1, 2);
  perm = setdigit(perm, 2, 1);
  //perm = setdigit(perm, 3, 3);
  hashdb patternset = hashdb(1<<3);
  patternset.add(perm);
  timestamp_t start_time = get_timestamp();
  cout<<"Avoid set: ";
  displayperm(perm);
  cout<<"Number of avoiders of size "<<permsize<<" is "<<countavoiders(patternset, maxpatternsize, permsize)<<endl;
  timestamp_t end_time = get_timestamp();
  cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  return 0;
}
					   
