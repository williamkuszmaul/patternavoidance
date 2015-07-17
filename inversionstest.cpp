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
#include "fastavoidance.h"
using namespace std;

// tests conjecture 20 of http://arxiv.org/abs/1111.5736 for permutations of size <= testsize avoiding patterns of size <= patternsize 

#define testsize 10
#define patternsize 7


int numinversions(uint64_t perm, int length) {
  int answer = 0;
  uint64_t seenletters = 0; // bit map of which letters we've seen so far
  for (int i = length - 1; i >= 0; i--) {
    if (getdigit(perm, i) != 0) answer += __builtin_popcountll(seenletters << (64 - getdigit(perm, i))); // note if digit = 0, the bit shifting is ill-defined
    seenletters = seenletters | (1 << getdigit(perm, i));
  }
  return answer;
}

int maxinversionsforsize(int size) {
  return size * (size - 1) / 2;
}

void fillinversions(hashdb &patternset, int avoidsize) {
  vector < vector < uint64_t > > avoidersvector;
  
  uint64_t numavoiders = buildavoiders(patternset, avoidsize, testsize, avoidersvector);

  int dim1 = testsize + 1, dim2 = maxinversionsforsize(testsize) + 1;
  int inversiontable[dim1][dim2]; // [i][j] stores num permutations of size i with inversion number j
  for (int i = 0; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      inversiontable[i][j] = 0;
    }
  }

  for (int permsize = 1; permsize < dim1; permsize++) {
    for (int permindex = 0; permindex < avoidersvector[permsize].size(); permindex++) {
      uint64_t avoider = avoidersvector[permsize][permindex];
      inversiontable[permsize][numinversions(avoider, permsize)]++;
    }
  }

  for (int i = 1; i < dim1; i++) {
    for (int j = 0; j < dim2; j++) {
      cout<<inversiontable[i][j]<<" ";
      assert(inversiontable[i-1][j] <= inversiontable[i][j]);
    }
    cout<<endl;
  }
  cout<<endl<<endl;
}

bool isid(uint64_t perm, int length) {
  for (int i = 1; i < length; i++) {
    if (getdigit(perm, i-1) > getdigit(perm, i)) return false;
  }
  return true;
}

void buildpermutations(uint64_t perm, int currentsize, int finalsize) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize);
    }
  } else {
    if (getdigit(perm, 0) == 0 && getdigit(perm, finalsize - 1) == finalsize - 1 && !isid(perm, finalsize)) {
      int avoidsize = finalsize;
      hashdb patternset = hashdb(1<<3);
      patternset.add(perm);
      cout<<"Inversion table for:"<<endl;
      displayperm(perm);
      cout<<endl;
      fillinversions(patternset, avoidsize);
      cout<<"--------------------------"<<endl;
    }
  }
}

int main() {
  uint64_t startperm = 0;
  buildpermutations(startperm, 1, patternsize);
  return 0;
}
