// g++ -O0 -std=c++11 -g countpatterns.cpp hashdb.h hashdb.cpp fastavoidance.h fastavoidance.cpp -o countpatterns
//
// scp countpatterns ec2:./
// ssh ec2
// /usr/bin/time ./countpatterns

// top gets things like gb usage. run in different shell
// can get big speedup by planning hash table size appropriately

//ssh ec2 to get into ec2
//scp countpatterns ec2:./
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

const int permsize = 5;
const int maxpatternsize = 4;


template <unsigned int MAXPATTERNSIZE> 
class myarray {
  public:
    unsigned short dataset[MAXPATTERNSIZE];
    myarray() {}
    myarray(unsigned short *array) {
    for (int i = 0 ; i < MAXPATTERNSIZE; i++) {
      dataset[i] = array[i];
    }
  }
};


inline void setPvals(unsigned long long perm, unsigned short* payload, unordered_map<unsigned long long, myarray<maxpatternsize + 1>> &Pmap) {
  myarray <maxpatternsize + 1> temp(payload);
  Pmap[perm] = temp;
}

inline unsigned short getPval(unsigned long long perm, int i, unordered_map<unsigned long long, myarray<maxpatternsize + 1>> &Pmap) {
  return Pmap[perm].dataset[i];
}

static void Pcount(uint64_t perm, int length, unordered_map<unsigned long long, myarray<maxpatternsize + 1>> &Pmap, const hashdb &patternset) { 
  uint64_t inverse = getinverse(perm, length);
  unsigned short Pvals[maxpatternsize + 2]; // vals range from [0...maxpatternsize+1]
  int Pvalpos = maxpatternsize + 1;
  while (Pvalpos > maxpatternsize || Pvalpos > length) {
    Pvals[Pvalpos] = 0;
    Pvalpos--;
  }
  if (Pvalpos == length) {
    if (patternset.contains(perm)) Pvals[Pvalpos] = 1;
    else Pvals[Pvalpos] = 0;
    Pvalpos--;
  }

  // Now Pvalpos < length and <= maxpatternsize. So the recurrence starts to apply

  unsigned short oldPvals[Pvalpos + 1];
  int oldPvalpos = 0;
  uint64_t currentperm = perm;
  assert(length > 1); // don't want to deal with length 1 case
  for (int i = length - 1; i >= 0 && i >= length - maxpatternsize - 1; i--) {
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
  setPvals(perm, Pvals, Pmap);
  displayperm(perm);
  cout<<"num hits: "<<Pvals[0]<<endl;
}


void buildpermutations(uint64_t perm, int currentsize, int finalsize, unordered_map<unsigned long long, myarray<maxpatternsize + 1>> &Pmap, hashdb &patternset) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize, Pmap, patternset);
    }
  } else {
    Pcount(perm, currentsize, Pmap, patternset);
  }
}

// LENGTH IS FUNNY SOME TIMES!
void createPmap(uint64_t finalsize, hashdb &patternset, unordered_map<unsigned long long, myarray<maxpatternsize + 1>> &Pmap, timestamp_t start_time) {
  unsigned short temp[4] = {0, 0, 0, 0};
  setPvals(0, temp, Pmap);
  uint64_t perm1 = setdigit(0L, 0, 1L);
  uint64_t perm2 = setdigit(0L, 1, 1L);
  for (int i = 2; i <= finalsize; i++) {
    buildpermutations(0L, 1, i, Pmap, patternset);
    //timestamp_t current_time = get_timestamp();
    //cout<< "Time elapsed to build perms of size "<<i<<" in seconds: "<<(current_time - start_time)/1000000.0L<<endl;
  }
  return;
}

int main() {
  assert(permsize <= 16);
  uint64_t perm = 0;
  perm = setdigit(perm, 0, 0);
  perm = setdigit(perm, 1, 1);
  perm = setdigit(perm, 2, 2);
  perm = setdigit(perm, 3, 3);
  hashdb patternset = hashdb(1<<3);
  patternset.add(perm);
 
  unordered_map<unsigned long long, myarray<maxpatternsize + 1>> Pmap; 
  timestamp_t start_time = get_timestamp();
  cout<<"Pattern set: ";
  displayperm(perm);
  createPmap(permsize, patternset, Pmap, start_time);
  timestamp_t end_time = get_timestamp();
  cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  return 0;
}
