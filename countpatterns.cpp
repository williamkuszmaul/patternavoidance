// Uses algorithm described paper (also on this github)
// We will use similar notation as in paper.pdf
// However, rather than P_i(perm) being the number of patterns using the first i letters in perm, it is the number o patterns using the i smallest-valued letters in perm (is equivalent from algorithmic perspective)

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
#include "perm.h"
using namespace std;

// We store P_0(perm), ..., P_{maxpatternsize}(perm) in a hash table containing arrays of shorts of size maxpatternsize + 1
inline void setPvals(unsigned long long perm, unsigned short* payload, hashmap &Phashmap) {
  Phashmap.add(perm, payload);
}

// Returns P_i(perm)
inline unsigned short getPval(unsigned long long perm, int i, hashmap &Phashmap) {
  return ((unsigned short*)Phashmap.getpayload(perm))[i];
}

// updates tally
void increasetally(vector < int > &tally, uint64_t val) {
  if (val + 1 > tally.size()) {
    uint64_t oldsize = tally.size();
    tally.resize(oldsize * 2);
    for (int i = oldsize; i < oldsize * 2; i++) tally[i] = 0;
  }
  tally[val]++;
}

// Inputs perm \in S_length, set of pattern patternset with longest pattern of size maxpatternsize, tally, completelist.
// Computes P_i(perm) for i from 0 , ... , maxpatternsize + 1; and if length < maxpermsize. stores all but final one
// Updates tally and completelist using P_0(perm)
static void Pcount(uint64_t perm, int length, int maxpatternsize, int maxpermsize, const hashdb &patternset, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist) {
  uint64_t inverse = getinverse(perm, length);
  unsigned short Pvals[maxpatternsize + 2]; // vals range from [0...maxpatternsize+1]
  int Pvalpos = maxpatternsize + 1;

  // Fill in trivially 0 Pvals
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
  unsigned short oldPvals[Pvalpos + 1]; // will store P_i(w\downarrow_{i+1}). We will then use these values for recurrence
  int oldPvalpos = 0;
  uint64_t currentperm = perm;
  assert(length > 1); // don't want to deal with length 1 case
  for (int i = length - 1; i >= 0 && i >= length - maxpatternsize - 1; i--) {
    if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
      currentperm = addpos(currentperm, getdigit(inverse, i + 1));
      currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
    }
    currentperm = killpos(currentperm, getdigit(inverse, i)); // now currentperm is perm, except with i removed and the resulting permutation normalized to be on values 0, ..., length - 1
    oldPvals[oldPvalpos] = getPval(currentperm, length - 1 - i, Phashmap);
    oldPvalpos++;
  }

  while (Pvalpos >= 0) {
    Pvals[Pvalpos] = Pvals[Pvalpos + 1] + oldPvals[Pvalpos]; // use recurrence to get actual Pvals
    Pvalpos--;
  }
  if (length != maxpermsize) {
    setPvals(perm, Pvals, Phashmap);
  }
  int index = permtonum(perm, length);
  completelist[length][index] = Pvals[0];
  increasetally(tally[length], Pvals[0]);
  //  displayperm(perm);
  //cout<<"num hits: "<<Pvals[0]<<endl;
}


// constructs the permutations of size finalsize. Then passes on paramaters to Pcount, in order to find number of patterns appearing in each permutation
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

// Fills in tally, complete list, and Pvals for all permutation in S_{<= finalsize}
// Note: patterns in patternset required to be in S_{>1}
void createPmap(uint64_t finalsize, hashdb &patternset, int maxpatternsize, timestamp_t start_time, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, bool verbose) {
  unsigned short temp[maxpatternsize + 1];
  for (int i = 0; i < maxpatternsize + 1; i++) temp[i] = 0;
  setPvals(0L, temp, Phashmap); // fill in Pvals to be 0 for 
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

// Example usage:
// string patternlist = "123 3124";
// int maxpermize = 10;
// vector < vector <int> > tally;
// vector < vector <int > > completelist;
// countpatterns(patternlist, maxpermsize, tally, completelist, false);
// Now tally[i][j] contains the number of permutation in S_i containing exactly j pattern-list-patterns. Range: 0 < i < maxpermsize + 1, 0 < j <= largest j such that tally[i][j] should exceed zero
// Now completelist[i][permtonum(perm)] is the number of pattern-set-hits appearing in perm, where perm \in S_i. Range 0 < i < maxsize + 1.
// Note: patterns in patternset required to be in S_{>1}
void countpatterns(string patternlist, int maxpermsize, vector < vector <int> > & tally, vector < vector < int > > &completelist, bool verbose) {
  assert(maxpermsize <= 16);

  // Start by resizing vectors appropriately
  tally.resize(maxpermsize + 1);
  completelist.resize(maxpermsize + 1);
  int fac = 1;
  for (int i = 1; i <= maxpermsize; i++) {
    fac *= i;
    completelist[i].resize(fac);
    tally[i].resize((1L << 10), 0);
  }
  // Note: tally's components are not appropriately sized at this point. This is done during the running of the code

  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize); // build pattern set

  unsigned long long reservedspace = 0;
  for (int i = 1; i <= maxpermsize - 1; i++) reservedspace += factorial(i);
  hashmap Phashmap(reservedspace * 3, sizeof(short)*(maxpatternsize + 1)); // initialize hash table of Pvals
  timestamp_t current_time = get_timestamp();
  createPmap(maxpermsize, patternset, maxpatternsize, current_time, Phashmap, tally, completelist, verbose);
  return;
}

// Example usage:
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
