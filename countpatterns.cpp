// g++ -O0 -std=c++11 -g countpatterns.cpp hashdb.h hashdb.cpp fastavoidance.h fastavoidance.cpp -o countpatterns

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
#include <unordered_map>
#include <queue>
#include "hashdb.h"
#include <sys/time.h>
#include "fastavoidance.h"
using namespace std;

#define testsize 11
#define patternsize 6

inline int setPval(unsigned long long perm, int i, int payload, unordered_map<unsigned long long, int> &Pmap) {
  std::pair<unsigned long long, int> pair ((perm<<4) + i, payload);
  Pmap.insert(pair);
}

inline int getPval(unsigned long long perm, int i, unordered_map<unsigned long long, int> &Pmap) {
  unordered_map<unsigned long long, int>::const_iterator got = Pmap.find((perm<<4) + i);
  assert(got != Pmap.end());
  return got->second;
}

static bool Pcount(uint64_t perm, int maxavoidsize, int length, unordered_map<unsigned long long, int> &Pmap, const hashdb &patternset) { 
  uint64_t inverse = getinverse(perm, length);
  int Pvals[maxavoidsize + 2]; // vals range from [0...maxavoidsize+1]
  int Pvalpos = maxavoidsize + 1;
  while (Pvalpos > maxavoidsize || Pvalpos > length) {
    Pvals[Pvalpos] = 0;
    Pvalpos--;
  }
  if (Pvalpos == length) {
    if (patternset.contains(perm)) Pvals[Pvalpos] = 1;
    else Pvals[Pvalpos] = 0;
    Pvalpos--;
  }

  // Now Pvalpos < length and <= maxavoidsize. So the recurrence starts to apply

  int oldPvals[Pvalpos + 1];
  int oldPvalpos = 0;
  uint64_t currentperm = perm;
  assert(length > 1); // don't want to deal with length 1 case
  for (int i = length - 1; i >= 0 && i >= length - maxavoidsize - 1; i--) {
    if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
      currentperm = addpos(currentperm, getdigit(inverse, i + 1));
      currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
    }
    currentperm = killpos(currentperm, getdigit(inverse, i));
    oldPvals[oldPvalpos] = getPval(currentperm, length - 1 - i, Pmap);
    oldPvalpos++;
  }

  while (Pvalpos >= 0) {
    Pvals[Pvalpos] = Pvals[Pvalpos + 1] + oldPvals[Pvalpos];
    Pvalpos--;
  }
  for (int i = 0; i <= maxavoidsize; i++) {
    setPval(perm, i, Pvals[i], Pmap);
  }
  //displayperm(perm);
  //cout<<"num hits: "<<Pvals[0]<<endl;
}


void buildpermutations(uint64_t perm, int currentsize, int finalsize, unordered_map<unsigned long long, int> &Pmap, int maxavoidsize, hashdb &patternset) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize, Pmap, maxavoidsize, patternset);
    }
  } else {
    Pcount(perm, maxavoidsize, currentsize, Pmap, patternset);
  }
}

void createPmap(uint64_t finalsize, int maxavoidsize, hashdb &patternset, unordered_map<unsigned long long, int> &Pmap, timestamp_t start_time) {
  setPval(0, 0, 0, Pmap);
  setPval(0, 1, 0, Pmap);
  uint64_t perm1 = setdigit(0L, 0, 1L);
  uint64_t perm2 = setdigit(0L, 1, 1L);
  for (int i = 2; i <= finalsize; i++) {
    buildpermutations(0L, 1, i, Pmap, maxavoidsize, patternset);
    timestamp_t current_time = get_timestamp();
    cout<< "Time elapsed to build perms of size "<<i<<" in seconds: "<<(current_time - start_time)/1000000.0L<<endl;
  }
  return;
}

int main() {
  int maxpatternsize = 3;
  int permsize = 12;
  assert(permsize <= 16);
  uint64_t perm = 0;
  perm = setdigit(perm, 0, 0);
  perm = setdigit(perm, 1, 2);
  perm = setdigit(perm, 2, 1);
  hashdb patternset = hashdb(1<<3);
  patternset.add(perm);
 
  unordered_map<unsigned long long, int> Pmap; 
  timestamp_t start_time = get_timestamp();
  cout<<"Pattern set: ";
  displayperm(perm);
  createPmap(permsize, maxpatternsize, patternset, Pmap, start_time);
  timestamp_t end_time = get_timestamp();
  cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  return 0;
}
