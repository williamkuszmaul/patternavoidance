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
#include <unordered_set>
#include <queue>
using namespace std;

#define ULL (uint64_t)
#define one ((uint64_t) 1)
#define zero ((uint64_t) 0)
#define bits4 ((uint64_t) 15)

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

bool isavoider(uint64_t perm, int maxavoidsize, int length, std::unordered_set<long long> &avoidset, std::unordered_set<long long> &patternset) { 
  uint64_t inverse = getinverse(perm, length);
  if (patternset.find(perm) != patternset.end()) { // if is in set of bad patterns
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
      if (avoidset.find(currentperm)  == patternset.end()) {
	return false; // found a subword not avoiding the patterns
      }
    }
  }
  return true;
}

int countavoiders(unordered_set<long long> patternset, int maxavoidsize, int maxsize) {
  int numavoiders = 0;
  
  std::unordered_set<long long> avoidset;
  uint64_t startperm = 0;
  avoidset.insert(0); // identity in S_1
  
  
  std::queue<long long> avoiderstoextend;
  avoiderstoextend.push(startperm);

  int currentlength = 1; // maintain as length of next permutation to be popped from avoiderstoextend
  int numleftcurrentlength = 1; // number of permutations left in avoiderstoextend until we have to increment currentlength
  int numnextlength = 0; // number of permutations of size currentlength + 1 in avoiderstoextend

  while (avoiderstoextend.empty() == false) {
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
	avoidset.insert(extendedperm);
      }
    }
  }
  return numavoiders;
}




int main() {
  uint64_t perm = 0;
  perm = setdigit(perm, 0, 0);
  perm = setdigit(perm, 1, 2);
  perm = setdigit(perm, 2, 1);
  perm = setdigit(perm, 3, 3);
  std::unordered_set<long long> patternset;
  patternset.insert(perm);
  cout<<countavoiders(patternset, 4, 11)<<endl;

  //for (int i = 0; i < 10; i++) perm = setdigit(perm, i, i);
  //isavoider(perm, 4, 6);
  
  
  // displayperm(perm);
  // displayperm(addpos(perm, 3));
  //displayperm(killpos(perm, 7));
}
