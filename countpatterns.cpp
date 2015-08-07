// g++ -O0 -std=c++11 -g countpatterns.cpp hashdb.h hashdb.cpp fastavoidance.h fastavoidance.cpp -o countpatterns
// Note: Valgrind falsely reports 72,704 bytes in use at exit, even when I return at beginning of main. This seems to be a documented valgrind bug.

// scp countpatterns ec2:./
// ssh ec2
// /usr/bin/time ./countpatterns

// top gets things like gb usage. run in different shell
// can get big speedup by planning hash table size appropriately

//ssh ec2 to get into ec2
//scp countpatterns ec2:./
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

unsigned long long permtonum(uint64_t perm, int length) {
  int answer = 0;
  unsigned long long fac = 1;
  uint64_t seenletters = 0; // bit map of which letters we've seen so far
  for (int i = length - 1; i >= 0; i--) {
    unsigned long long facdig = 0;
    int digit = getdigit(perm, i);
    if (digit != 0) {
      unsigned int temp = seenletters << (32 - digit);
      facdig = __builtin_popcount(temp);
    }
    answer += facdig * fac;
    fac *= (length - i);
    seenletters = seenletters | (1 << digit);
  }
  return answer;
}



inline void setPvals(unsigned long long perm, unsigned short* payload, hashmap &Phashmap) {
  Phashmap.add(perm, payload);
}

inline unsigned short getPval(unsigned long long perm, int i, hashmap &Phashmap) {
  return ((unsigned short*)Phashmap.getpayload(perm))[i];
}


void increasetally(vector < int > &tally, uint64_t val) {
  if (val + 1 > tally.size()) {
    uint64_t oldsize = tally.size();
    tally.resize(oldsize * 2);
    for (int i = oldsize; i < oldsize * 2; i++) tally[i] = 0;
  }
  tally[val]++;
}
				       
static void Pcount(uint64_t perm, int length, int maxpatternsize, int maxpermsize, const hashdb &patternset, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist) {
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
    oldPvals[oldPvalpos] = getPval(currentperm, length - 1 - i, Phashmap);
    //assert(getPval(currentperm, length - 1 - i, Pmap) ==
    //   getPval2(currentperm, length - 1- i, Phashmap));
    oldPvalpos++;
  }

  while (Pvalpos >= 0) {
    Pvals[Pvalpos] = Pvals[Pvalpos + 1] + oldPvals[Pvalpos];
    Pvalpos--;
  }
    if (length != maxpermsize) {
      //setPvals(perm, Pvals, Pmap);
      setPvals(perm, Pvals, Phashmap);
    }
    int index = permtonum(perm, length);
    completelist[length][index] = Pvals[0];
    increasetally(tally[length], Pvals[0]);
  //  displayperm(perm);
  //cout<<"num hits: "<<Pvals[0]<<endl;
}


void buildpermutations(uint64_t perm, int currentsize, int finalsize, int maxpatternsize, int maxpermsize, hashdb &patternset, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize, maxpatternsize,  maxpermsize, patternset, Phashmap, tally, completelist);
    }
  } else {
    Pcount(perm, currentsize, maxpatternsize, maxpermsize, patternset, Phashmap, tally, completelist);
  }
}


void createPmap(uint64_t finalsize, hashdb &patternset, int maxpatternsize, timestamp_t start_time, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, bool verbose) {
  unsigned short temp[4] = {0, 0, 0, 0};
  setPvals(0, temp, Phashmap);
  uint64_t perm1 = setdigit(0L, 0, 1L);
  uint64_t perm2 = setdigit(0L, 1, 1L);
  for (int i = 2; i <= finalsize; i++) {
    buildpermutations(0L, 1, i, maxpatternsize, finalsize, patternset, Phashmap, tally, completelist);
    timestamp_t current_time = get_timestamp();
    if (verbose) cout<< "Time elapsed to build perms of size "<<i<<" in seconds: "<<(current_time - start_time)/1000000.0L<<endl;
  }
  return;
}

unsigned long long factorial(long long num) {
  int answer = 1;
  while (num > 1) {
    answer *= num;
    num--;
  }
  return answer;
}

void countpatterns(string patternlist, int maxpermsize, vector < vector <int> > & tally, vector < vector < int > > &completelist, bool verbose) {
  assert(maxpermsize <= 16);
  
  tally.resize(maxpermsize + 1);
  completelist.resize(maxpermsize + 1);
  int fac = 1;
  for (int i = 1; i <= maxpermsize; i++) {
    fac *= i;
    completelist[i].resize(fac);
    tally[i].resize((1L << 10), 0);
  }
  // Note tally's components are not appropriately sized at this point

  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);

  unsigned long long reservedspace = 0;
  for (int i = 1; i <= maxpermsize - 1; i++) reservedspace += factorial(i);
  hashmap Phashmap(reservedspace * 3, sizeof(short)*(maxpatternsize + 1));
  timestamp_t current_time = get_timestamp();
  createPmap(maxpermsize, patternset, maxpatternsize, current_time, Phashmap, tally, completelist, verbose);
  return;
}

// int main() {
//   vector < vector < int > > tally;
//   vector < vector < int > > completelist;
//   countpatterns("1234 123", 8, tally, completelist);
//   for (int i = 1; i <= 8; i ++) cout<<completelist[i][0]<<endl;
//   cout<<"-----------------"<<endl;
//   for (int i = 8; i <= 8; i++) {
//     for (int j = 0; j < tally[i].size(); j++) {
//       cout<<j<<" "<<tally[i][j]<<endl;
//     }
//   }
//   return 0;
// }
