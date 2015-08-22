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

#define USEADDFACTOR 1
//A#define USESECONDADDFACTOR 1 // Requires USEADDFACTOR on

// Input: perm, inverse (which needs to be right in pos length - index - 1), perm length, index, complement of normalization of index - 1 largest letters in perm, a bitmap which should start off at zero for index = 0
// Output: bitmap is updated. answer is updated to be complement of normalization of index largest letters in perm
static void extendnormalizetop(uint64_t perm, uint64_t inverse, int length, int index, uint64_t &answer, uint32_t & seenpos) {
  int i = length - index - 1;
  int oldpos = getdigit(inverse, i);
  int newpos = 0;
  if (oldpos != 0){
    uint32_t temp = seenpos << (32 - oldpos); // Note: shifting by 32 is ill-defined, which is why we explicitly eliminate digit = 0 case.
    newpos = __builtin_popcount(temp);
  }
  answer = setdigit(addpos(answer, newpos), newpos, index);
  seenpos = seenpos | (1 << oldpos);
}


static void addprefixeshelper(uint64_t perm, int length, hashdb &table) {
  uint64_t entry = 0;
  uint64_t inverse = getinverse(perm, length);
  uint32_t seenpos = 0; // bit map of which letters we've seen so far
  for (int i = 0; i < length; i++) {
    extendnormalizetop(perm, inverse, length, i, entry, seenpos);
    //cout<<entry<<endl;
    if (!table.contains(entry)) table.add(entry);
  }
}

// build prefix table
static void addprefixes(const hashdb &permset, hashdb &table) {
  vector <unsigned long long> patterns;
  permset.getvals(patterns);
  for (int i = 0; i < patterns.size(); i++) {
    uint64_t perm = patterns[i];
    int length = getmaxdigit(perm) + 1;
    displayperm(perm,length);
    addprefixeshelper(perm, length, table);
  }
}


void checkpatterns(uint64_t perm, uint64_t inverse, uint64_t currentpatterncomplement, int currentpatternlength, int largestletterused, int numlettersleft, uint32_t seenpos, const hashdb &patternset, const hashdb &prefixmap, int &count) {
  if (currentpatterncomplement != 0 && !prefixmap.contains(currentpatterncomplement)) return;
  if (currentpatterncomplement != 0 && patternset.contains(currentpatterncomplement)) count++;
  if (numlettersleft == 0) return; // make sure this if statement comes AFTER checking for patternset
  for (int i = largestletterused - 1; i >= 0; i--) {
    if (USEADDFACTOR && currentpatternlength == 0 && i < largestletterused - 1) return; // because of our use of addfactor
    //if (USESECONDADDFACTOR && currentpatternlength == 1 && i < largestletterused - 2) return; 
    int oldpos = getdigit(inverse, i);
    int newpos = 0;
    if (oldpos != 0){  
      uint32_t temp = seenpos << (32 - oldpos); // Note: shifting by 32 is ill-defined, which is why we explicitly eliminate digit = 0 case.
      newpos = __builtin_popcount(temp);
    }
    //cout<<"Here again with i selected at "<<i<<endl;
    checkpatterns(perm, inverse, setdigit(addpos(currentpatterncomplement, newpos), newpos, currentpatternlength), currentpatternlength + 1, i, numlettersleft - 1, seenpos | (1 << oldpos), patternset, prefixmap, count);
  }
  return;
}

// constructs the permutations of size finalsize. Then passes on paramaters to Pcount, in order to find number of patterns appearing in each permutation
void buildpermutations_brute(uint64_t perm, uint64_t inverse, int currentsize, int finalsize, int maxpatternsize, int maxpermsize, hashdb &patterncomplements, hashdb &prefixmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, int addfactor) {
      //    uint64_t inverse = getinverse(perm, currentsize);
  int count = 0;
  checkpatterns(perm, inverse, 0, 0, currentsize, maxpatternsize, 0, patterncomplements, prefixmap, count);
  if (USEADDFACTOR) count += addfactor;
  tally[currentsize][count]++;
  completelist[currentsize][permtonum(perm, currentsize)] = count;
  
  uint64_t newinverse = setdigit(inverse, currentsize, currentsize); // inverse of the extended permutation
  if (currentsize < finalsize) {
    for (int i = currentsize; i >= 0; i--) {
      if (i < currentsize) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * currentsize));
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      assert(getinverse(extendedperm, currentsize + 1) == newinverse);
      buildpermutations_brute(extendedperm, newinverse, currentsize + 1, finalsize, maxpatternsize,  maxpermsize, patterncomplements, prefixmap, tally, completelist, count);
    }
  }
}

void start_brute(int maxpatternsize, int maxpermsize, hashdb &patternset, vector < vector < int > > &tally, vector < vector < int > > &completelist) {
  cout<<"Using Brute Force Algorithm"<<endl;

  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap);
  
  vector <unsigned long long> patterns;
  hashdb patterncomplements(1<<3);
  patternset.getvals(patterns);
  for (int i = 0; i < patterns.size(); i++) {
    uint64_t perm = patterns[i];
    int length = getmaxdigit(perm) + 1;
    patterncomplements.add(getcomplement(perm, length));
  }
  cout<<"max size: "<<maxpermsize<<endl;
  for (int i = 2; i <= maxpermsize; i++) {
    buildpermutations_brute(0L, 0L, 1, i, maxpatternsize, maxpermsize, patterncomplements, prefixmap, tally, completelist, 0);
    timestamp_t current_time = get_timestamp();
    //    if (verbose) cout<< "Time elapsed to build perms of size "<<i<<" in seconds: "<<(current_time - start_time)/1000000.0L<<endl;
  }
}

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
  start_brute(maxpatternsize, maxpermsize, patternset, tally, completelist);
  //createPmap(maxpermsize, patternset, maxpatternsize, current_time, Phashmap, tally, completelist, verbose);
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
