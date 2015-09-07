
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

// Conventions for comments: We use all the conventions from our paper in our comments. Thus the comments disagree with the code in that they do not zero-index the values or positions of a permutation

// Turn these on to turn brute force algorithm into a memory efficient hybrid between the brute-force algorithm and the asymptotically good algorithm
#define USEADDFACTOR 1
#define USESECONDADDFACTOR 1 // forces useaddfactor on

#define USEOLDP0 1
#define USEOLDP1 1
#define USEPREFIXMAP 1 // FOR NON-BRUTE FORCE VERSION
#define USEBRUTE 0


unsigned long long factorial(long long num) {
  int answer = 1;
  while (num > 1) {
    answer *= num;
    num--;
  }
  return answer;
}

// updates tally
void increasetally(vector < int > &tally, uint64_t val) {
  if (val + 1 > tally.size()) {
    uint64_t oldsize = tally.size();
    uint64_t size = oldsize;
    while (size < val + 1) {
      size *= 2;
      tally.resize(size); // resize if needed
    }
    for (int i = oldsize; i < size; i++) tally[i] = 0;
  }
  tally[val]++;
}

// Recursively checks whether perm contains a pattern from patternset.
// Requires that all patterns in patternset are length currentpatternlength + numlettersleft
// currentpatterncomplement stores complement of normalization of permutation subsequence already examined
// largestletterused tells us the value of the largest letter icnluded so far in the permutation-subsequence being examined
// seenpos is a bitmap used to efficiently update currentpatterncomplement as we make the subsequence-being looked at one longer
// prefixmap contains the complements of each prefix of each pattern in \Pi.
// Note that prefixmap contains complements of normalized prefixes of patterns
// If USEADDFACTOR, computes P_1(perm) instead of P(perm). If USESECONDADDFACTOR, computes P_2(perm) instead
static void checkpatterns(uint64_t perm, uint64_t inverse, uint64_t currentpatterncomplement, int currentpatternlength, int largestletterused, int numlettersleft, uint32_t seenpos, const hashdb &prefixmap, int &count) {
  
  if (currentpatterncomplement != 0 && !prefixmap.contains(currentpatterncomplement)) return;
  if (numlettersleft == 0) { // Assumes all patterns are same size --> this is only cae where prefix is a pattern
    count++;
    return; // also assumes all patterns same size
  }
  if (numlettersleft == 0) return;
  for (int i = largestletterused - 1; i >= numlettersleft - 1; i--) { // for i >= numlettersleft requirement, assumes all patterns are same size 
    if ((USEADDFACTOR || USESECONDADDFACTOR) && currentpatternlength == 0 && i < largestletterused - 1) return; // because of our use of addfactor
    if (USESECONDADDFACTOR && currentpatternlength == 1 && i < largestletterused - 1) return; // note this means i < maxpatternlength - 2, because if-statement above forces largestletterused to decrease by 1 to maxpatternlength - 1 in the previous level of recursion

    // Similarly to as in extendnormalizetop (defined in perm.cpp), we will build the complement of the normalization of the new permutation-subsequence (with the new letter added to it)
    int oldpos = getdigit(inverse, i);
    int newpos = 0;
    if (oldpos != 0){  
      uint32_t temp = seenpos << (32 - oldpos); //Note: shifting by 32 is ill-defined, which is why we explicitly eliminate digit = 0 case.
      newpos = __builtin_popcount(temp);
    }
    uint64_t newpattern = setdigit(addpos(currentpatterncomplement, newpos), newpos, currentpatternlength);
    // Recurse to make sequence longer until we eventually either do or don't get a pattern:
    checkpatterns(perm, inverse, newpattern, currentpatternlength + 1, i, numlettersleft - 1, seenpos | (1 << oldpos), prefixmap, count);
  }
  return;
}

// constructs the permutations of size finalsize. Then passes on paramaters to checkpatterns, in order to find number of patterns appearing in each permutation
void buildpermutations_brute(uint64_t perm, uint64_t inverse, int currentsize, int finalsize, int maxpatternsize, hashdb &prefixmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, int addfactor, bool justcount) {
  int count = 0;
  checkpatterns(perm, inverse, 0, 0, currentsize, maxpatternsize, 0, prefixmap, count);
  if (USEADDFACTOR) count += addfactor; // add factor is P_0(perm\downarrow_1)
  increasetally(tally[currentsize], count);
  if (!justcount) completelist[currentsize][permtonum(perm, currentsize)] = count;
  
  uint64_t newinverse = setdigit(inverse, currentsize, currentsize); // inverse of perm\uparrow^{currentsize+1}
  if (currentsize < finalsize) {
    for (int i = currentsize; i >= 0; i--) {
      // need to increment newinverse[perm[i]], decrement newinverse[length]
      if (i < currentsize) newinverse = decrementdigit(incrementdigit(newinverse, getdigit(perm, i)), currentsize);
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize); // extendedperm = perm\uparrow^{i + 1}
      buildpermutations_brute(extendedperm, newinverse, currentsize + 1, finalsize, maxpatternsize, prefixmap, tally, completelist, count, justcount);
    }
  }
}

// same as buildpermutations_brute but with USESECONDADDFACTOR installed
void buildpermutations_brute_usingbothaddfactors(uint64_t perm, uint64_t inverse, int currentsize, int finalsize, int maxpatternsize, hashdb &prefixmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, int prevcount, int* prevextensions, bool justcount) {
  int actualcounts[currentsize + 1]; // Each P_0(perm\uparrow^{i + 1})
  int newextensions[currentsize + 1]; // Each P_1(perm\uparrow^{i + 1}) 
  uint64_t newinverse = setdigit(inverse, currentsize, currentsize); // inverse of the perm\uparrow^{currentsize + 1}
  if (currentsize < finalsize) {
    for (int i = currentsize; i >= 0; i--) {
      if (i < currentsize) newinverse = decrementdigit(incrementdigit(newinverse, getdigit(perm, i)), currentsize);
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize); // is perm\uparrow^{i + 1}
      int countusingfinaltwo = 0; // will be P_2(extendedperm)
      checkpatterns(extendedperm, newinverse, 0, 0, currentsize + 1, maxpatternsize, 0, prefixmap, countusingfinaltwo); // fill in countusingfinaltwo
      int indexignoringsecondtofinal = i;
      if (getdigit(newinverse, currentsize) > getdigit(newinverse, currentsize - 1)) indexignoringsecondtofinal--; 
      int actualcount = countusingfinaltwo + prevcount + prevextensions[indexignoringsecondtofinal]; // P_0(extendedperm) = P_1(extendedperm) + P_0(extendedperm\downarrow_1) = (P_2(extendedperm) + P_1(extendedperm\downarrow_2)) + P_0(extendedperm\downarrow_1)
      increasetally(tally[currentsize + 1], actualcount);
      if (!justcount) completelist[currentsize + 1][permtonum(extendedperm, currentsize + 1)] = actualcount;
      newextensions[i] = countusingfinaltwo + prevextensions[indexignoringsecondtofinal]; // P_1(extendedperm) = P_2(extendedperm) + P_1(extendedperm\downarrow_2)
      actualcounts[i] = actualcount;
    }  
  }
  
  newinverse = setdigit(inverse, currentsize, currentsize); // inverse of the extended permutation
  if (currentsize < finalsize - 1) {
    for (int i = currentsize; i >= 0; i--) {
      // need to increment newinverse[perm[i]], decrement newinverse[length]
      if (i < currentsize) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * currentsize));
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations_brute_usingbothaddfactors(extendedperm, newinverse, currentsize + 1, finalsize, maxpatternsize, prefixmap, tally, completelist, actualcounts[i], newextensions, justcount);
    }
  }
}

void start_brute(int maxpatternsize, int maxpermsize, hashdb &patternset, vector < vector < int > > &tally, vector < vector < int > > &completelist, bool justcount) {
  cout<<"Using Brute Force Algorithm"<<endl;

  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap);

  
  vector <unsigned long long> patterns;
  //hashdb patterncomplements(1<<3);
  patternset.getvals(patterns);
  // for (int i = 0; i < patterns.size(); i++) {
  //   uint64_t perm = patterns[i];
  //   int length = getmaxdigit(perm) + 1;
  //   patterncomplements.add(getcomplement(perm, length));
  // }
  timestamp_t start_time = get_timestamp();
  cout<<"max size: "<<maxpermsize<<endl;
  if (!USESECONDADDFACTOR) {
    buildpermutations_brute(0L, 0L, 1, maxpermsize, maxpatternsize, prefixmap, tally, completelist, 0, justcount);
  } else {
    int temparray[2] = {0, 0};
    buildpermutations_brute_usingbothaddfactors(0L, 0L, 1, maxpermsize, maxpatternsize, prefixmap, tally, completelist, 0, temparray, justcount);
  }
  timestamp_t current_time = get_timestamp();
  cout<< "Time elapsed to build perms of size "<<maxpermsize<<" in seconds: "<<(current_time - start_time)/1000000.0L<<endl;
}

// We store P_0(perm), ..., P_{maxpatternsize}(perm) in a hash table containing arrays of shorts of size maxpatternsize + 1
inline void setPvals(unsigned long long perm, unsigned short* payload, hashmap &Phashmap) {
  Phashmap.add(perm, payload);
}

// Returns P_i(perm)
inline unsigned short getPval(unsigned long long perm, int i, hashmap &Phashmap) {
  //assert(Phashmap.getpayload(perm) != NULL);
  return ((unsigned short*)Phashmap.getpayload(perm))[i];
}


// Inputs perm \in S_length, set of pattern patternset with longest pattern of size maxpatternsize, tally, completelist.
// Computes P_i(perm) for i from 0 , ... , maxpatternsize + 1; and if length < maxpermsize. stores all but final one
// Updates tally and completelist using P_0(perm). Only updates completelist of justcount is false
// Phashmap contains all Pvals previously computed (must include everything for S_{length - 1}.
// patternset contains the complements of the normalizations of the prefixes of patterns
static void Pcount(uint64_t perm, int length, int maxpatternsize, int maxpermsize, const hashdb &patternset, const hashdb &prefixmap, hashmap &Phashmap, uint64_t prevP0, uint64_t prevP1,  vector < vector < int > > &tally, vector < vector < int > > &completelist, bool justcount) {
  uint64_t inverse = getinverse(perm, length);
  unsigned short Pvals[maxpatternsize + 2]; // indices range from [0...maxpatternsize+1]. Assumes no Pvals as short as large as 63504. Will not work, for example, for counting id_10 patterns in id_20.
  int Pvalpos = maxpatternsize + 1;

  // Don't have to worry about Pvalpos[i] if i > length
  Pvalpos = min(Pvalpos,  length);
  Pvals[Pvalpos] = 0;
  if (Pvalpos == length) {
    if (patternset.contains(perm)) Pvals[Pvalpos] = 1;
  }
  Pvalpos--;

  // Now Pvalpos < length and <= maxpatternsize. So the recurrence starts to apply
  unsigned short oldPvals[Pvalpos + 1]; // will store P_i(perm\downarrow_i). We will then use these values for recurrence
  int oldPvalpos = 0;
  uint64_t currentperm = perm; // will be updated to be each perm \downarrow_i
  assert(length > 1); // don't want to deal with length 1 case
  uint32_t seenpos = 0;
  uint64_t prefixentry = 0;

  for (int i = length - 1; i >= 0 && i >= length - maxpatternsize - 1; i--) {
    if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
      currentperm = addpos(currentperm, getdigit(inverse, i + 1));
      currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
    }
    currentperm = killpos(currentperm, getdigit(inverse, i)); // now currentperm is perm, except with i removed and the resulting permutation normalized to be on values 0, ..., length - 1
    if (USEOLDP0 && oldPvalpos == 0) {
      oldPvals[0] = prevP0;
    } else {
      if (USEOLDP1 && oldPvalpos == 1) {
	oldPvals[1] = prevP1;
      } else {
	oldPvals[oldPvalpos] = getPval(currentperm, length - 1 - i, Phashmap);
      }
    }
    oldPvalpos++;
    if (USEPREFIXMAP) {
      extendnormalizetop(perm, inverse, length, length - 1 - i, prefixentry, seenpos); // build the complement fo the normalized (length - i)-prefix of the permutation.
      if (!prefixmap.contains(prefixentry)) { // If this prefix is not a valid pattern prefix
	Pvalpos = oldPvalpos - 1; 
	Pvals[oldPvalpos] = 0; // we know that Pvals[oldPvalpos] is 0
	// We do not have to both filling in Pvals[i] for i > oldPvalpos, because there will never be a permutation $w$ such that $w \downarrow_i = perm$, i > oldPvalpos, and the (i-1)-prefix o $w$ forms a valid pattern prefix.
	// for (int j = oldPvalpos; j <= Pvalpos; j++) { // This is what we would have done instead of previous two lines if we did want to fill in all the Pvals
	//   oldPvals[j] = 0; 
	// }
	break;
      }
    }
  }

  // Now we start the recursion, starting with P_{Pvalpos})(perm)
  while (Pvalpos >= 0) {
    Pvals[Pvalpos] = Pvals[Pvalpos + 1] + oldPvals[Pvalpos]; // use recurrence to get actual Pvals
    Pvalpos--;
  }
  if (length != maxpermsize) {
    setPvals(perm, Pvals, Phashmap);
  }
  if (!justcount) {
    int index = permtonum(perm, length); // Note: permtonum is O(n) time
    completelist[length][index] = Pvals[0];
  }
  increasetally(tally[length], Pvals[0]);
}


// constructs the permutations of size finalsize. Then passes on paramaters to Pcount, in order to find number of patterns appearing in each permutation
// Requires that it has already been run for each smaller currentsize >= 1
// Note first call to this function should be with currentsize = 1
void buildpermutations(uint64_t perm, int currentsize, int finalsize, int maxpatternsize, int maxpermsize, hashdb &patternset, hashdb &prefixmap, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, uint64_t *cachedP1s, bool justcount) {
  if (USEOLDP1 && currentsize == finalsize - 2) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      cachedP1s[i] = getPval(extendedperm, 1, Phashmap); // store these Pvals now, and they will be used many times two levels down from now
    }
  }
  bool nseen = false; // whether or not we've seen currentsize in perm yet
  uint64_t currentP0 = 0; // only will be filled in if currentsize == finalsize - 1
  if (currentsize == finalsize - 1) currentP0 = getPval(perm, 0, Phashmap);
  for (int i = 0; i < currentsize + 1; i++) {
    if (i > 0 && getdigit(perm, i - 1) == currentsize - 1) nseen = true;
    uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
    if (currentsize < finalsize - 1) {
      buildpermutations(extendedperm, currentsize + 1, finalsize, maxpatternsize,  maxpermsize, patternset, prefixmap, Phashmap, tally, completelist, cachedP1s, justcount);
    } else {
      uint64_t prevP1 = cachedP1s[i];
      if (nseen) prevP1 = cachedP1s[i - 1];
      Pcount(extendedperm, currentsize + 1, maxpatternsize, maxpermsize, patternset, prefixmap, Phashmap, currentP0, prevP1, tally, completelist, justcount);
    }
  }
}

// Fills in tally, complete list, and Pvals for all permutation in S_{<= finalsize} (unles justcount, inwhich case does not fill in completelist)
// Note: patterns in patternset required to be in S_{>1}
void createPmap(uint64_t finalsize, hashdb &patternset, int maxpatternsize, timestamp_t start_time, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, bool verbose, bool justcount) {
  unsigned short temp[maxpatternsize + 1];
  for (int i = 0; i < maxpatternsize + 1; i++) temp[i] = 0;
  setPvals(0L, temp, Phashmap); // fill in Pvals to be 0 for

  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap);

  for (int i = 2; i <= finalsize; i++) {
    uint64_t cachedP1s[finalsize];
    for (int x = 0; x < finalsize; x++) cachedP1s[x] = 0; // initialize to zero used in first call to buildpermutations when i = 2
    buildpermutations(0L, 1, i, maxpatternsize, finalsize, patternset, prefixmap, Phashmap, tally, completelist, cachedP1s, justcount);
    timestamp_t current_time = get_timestamp();
    if (verbose) cout<< "Time elapsed to build perms of size "<<i<<" in seconds: "<<(current_time - start_time)/1000000.0L<<endl;
  }
  return;
}


double run_interior_experiment2(string patternlist, int maxpermsize) {
  timestamp_t start_time = get_timestamp();
  vector < vector <int> > tally;
  vector < vector < int > > completelist;
  bool verbose = false;
  assert(maxpermsize <= 16);

  // Start by resizing vectors appropriately
  tally.resize(maxpermsize + 1);
  completelist.resize(maxpermsize + 1);
  int fac = 1;
  for (int i = 1; i <= maxpermsize; i++) {
    fac *= i;
    //completelist[i].resize(fac);
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
  if (USEBRUTE)  start_brute(maxpatternsize, maxpermsize, patternset, tally, completelist, true);
  else createPmap(maxpermsize, patternset, maxpatternsize, current_time, Phashmap, tally, completelist, verbose, true);  
  timestamp_t end_time = get_timestamp();
  return (end_time - start_time)/1000000.0L;
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
void countpatterns(string patternlist, int maxpermsize, vector < vector <int> > & tally, vector < vector < int > > &completelist, bool verbose, bool justcount) {
  timestamp_t start_time = get_timestamp();
  assert(maxpermsize <= 16);

  // Start by resizing vectors appropriately
  tally.resize(maxpermsize + 1);
  completelist.resize(maxpermsize + 1);
  int fac = 1;
  for (int i = 1; i <= maxpermsize; i++) {
    fac *= i;
    //completelist[i].resize(fac); // not currently used
    tally[i].resize((1L << 10), 0);
  }
  // Note: tally's components are not appropriately sized at this point. This is done during the running of the code
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);

  makepatterns(patternlist, patternset, maxpatternsize); // build pattern set
  unsigned long long reservedspace = 0;
  for (int i = 1; i <= maxpermsize - 1; i++) reservedspace += factorial(i);
  hashmap Phashmap(reservedspace * 3, sizeof(short)*(maxpatternsize + 1)); // initialize hash table of Pvals // should be automatically optimized out in brute-force version
  timestamp_t current_time = get_timestamp();
  if (!justcount) {
    completelist.resize(maxpermsize + 1);
    for (uint64_t i = 0; i < maxpermsize + 1; i++) {
      completelist[i].resize(factorial(i));
    }
  }

  if (USEBRUTE)  start_brute(maxpatternsize, maxpermsize, patternset, tally, completelist, justcount);
  else createPmap(maxpermsize, patternset, maxpatternsize, current_time, Phashmap, tally, completelist, verbose, justcount);
  timestamp_t end_time = get_timestamp();
  if (verbose) cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
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
