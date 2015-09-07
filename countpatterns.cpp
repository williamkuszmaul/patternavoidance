
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

// Turn these on to turn brute force algorithm into a memory efficient hybrid betweenthe brute-force algorithm and the asymptotically good algorithm
#define USEADDFACTOR 1
#define USESECONDADDFACTOR 1 // forces useaddfactor on
// Input: perm, inverse (which needs to be right in pos length - index - 1), perm length, index, complement of normalization of index - 1 largest letters in perm, a bitmap which should start off at zero for index = 0
// Output: bitmap is updated. answer is updated to be complement of normalization of index largest letters in perm

#define USEOLDP0 1
#define USEOLDP1 1
#define USEPREFIXMAP 1 // FOR NON-BRUTE FORCE VERSION
#define SPECIALTEST 0 // for if single identity pattern case
#define USEBRUTE 1

static uint64_t stat1 = 0;
static uint64_t stat2 = 0;

// updates tally
void increasetally(vector < int > &tally, uint64_t val) {
  if (val + 1 > tally.size()) {
    uint64_t oldsize = tally.size();
    uint64_t size = oldsize;
    while (size < val + 1) {
      size *= 2;
      tally.resize(size);
    }
    for (int i = oldsize; i < size; i++) tally[i] = 0;
  }
  tally[val]++;
}


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
    //displayperm(perm,length);
    addprefixeshelper(perm, length, table);
  }
}


uint64_t lcount = 0;

void checkpatterns(uint64_t perm, uint64_t inverse, uint64_t currentpatterncomplement, int currentpatternlength, int largestletterused, int numlettersleft, uint32_t seenpos, const hashdb &patternset, const hashdb &prefixmap, int &count) {
  //uint64_t tempperm = stringtoperm("1234");
  
  if (SPECIALTEST) {
    if (numlettersleft == 0) {
      count++;
    }
  } else {
    if (currentpatterncomplement != 0 && !prefixmap.contains(currentpatterncomplement)) return;
    if (numlettersleft == 0) { // ASSUMES ALL PATTERNS SAME SIZE
      count++;
      return;
    }
    // for if not all patterns same size:
    // if (currentpatterncomplement != 0 && patternset.contains(currentpatterncomplement)) {
    //   count++;
    //   //if (perm == tempperm) cout<<"Here! top used letter is "<<largestletterused<<endl;
    // }
  }
  if (numlettersleft == 0) return; // make sure this if statement comes AFTER checking for patternset
  for (int i = largestletterused - 1; i >= numlettersleft - 1; i--) {
    if ((USEADDFACTOR || USESECONDADDFACTOR) && currentpatternlength == 0 && i < largestletterused - 1) return; // because of our use of addfactor
    if (USESECONDADDFACTOR && currentpatternlength == 1 && i < largestletterused - 1) return; // note this means i < maxpatternlength - 2, because if-statement above forces largestletterused to decrease by 1 to maxpatternlength - 1 in the previous level of recursion 
    int oldpos = getdigit(inverse, i);
    int newpos = 0;
    if (oldpos != 0){  
      uint32_t temp = seenpos << (32 - oldpos); //Note: shifting by 32 is ill-defined, which is why we explicitly eliminate digit = 0 case.
      newpos = __builtin_popcount(temp);
    }
    //cout<<"Here again with i selected at "<<i<<endl;
    //if (perm == tempperm) cout<<(currentpatternlength + 1)<<" to right-most letter is in position "<<i<<endl;
    uint64_t newpattern = setdigit(addpos(currentpatterncomplement, newpos), newpos, currentpatternlength);
    if (false && SPECIALTEST) { // If I don't run this, then we get almost exactly same running time for the sophisticated brute force and nonbrute force. So close, that it may not be coincidental since they may be doing almost exactly the same thing in this case; probably just coincidence though since getting rid fo division operation will probably speed up the non-brute-force version a bit.. Raises a natural question: Can this be extendd to small set of patterns case without having a few load operations total everything so that we may as well be using hash tble again. O
      bool cont = prefixmap.simulatelookup(newpattern);
      bool cont3 = prefixmap.simulatelookup(newpattern + 1);
      bool cont2 = false;
      if (cont) cont2 = patternset.simulatelookup(newpattern);
      static int count_cont = 0;
      count_cont += cont + cont2 + cont3;
      if (count_cont % 1024*1024*1024 == 9999999) cout << "whahoo" << endl;
      lcount++;
    }
    if (!SPECIALTEST || currentpatternlength == 0 || getdigit(inverse, i) < getdigit(inverse, largestletterused)) checkpatterns(perm, inverse, newpattern, currentpatternlength + 1, i, numlettersleft - 1, seenpos | (1 << oldpos), patternset, prefixmap, count);
  }
  return;
}

uint64_t tempcount = 0; // just for code testing ue
// constructs the permutations of size finalsize. Then passes on paramaters to Pcount, in order to find number of patterns appearing in each permutation
void buildpermutations_brute(uint64_t perm, uint64_t inverse, int currentsize, int finalsize, int maxpatternsize, int maxpermsize, hashdb &patterncomplements, hashdb &prefixmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, int addfactor) {
  int count = 0;
  checkpatterns(perm, inverse, 0, 0, currentsize, maxpatternsize, 0, patterncomplements, prefixmap, count);
  //tempcount++;
  //if (tempcount == 1000000L) cout<<tempcount<<endl;;
  if (USEADDFACTOR) count += addfactor;
  increasetally(tally[currentsize], count);
  //completelist[currentsize][permtonum(perm, currentsize)] = count;
  
  uint64_t newinverse = setdigit(inverse, currentsize, currentsize); // inverse of the extended permutation
  if (currentsize < finalsize) {
    for (int i = currentsize; i >= 0; i--) {
      if (i < currentsize) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * currentsize));
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations_brute(extendedperm, newinverse, currentsize + 1, finalsize, maxpatternsize,  maxpermsize, patterncomplements, prefixmap, tally, completelist, count);
    }
  }
}

// same as buildpermutations_brute but with USESECONDADDFACTOR installed
void buildpermutations_brute_usingbothaddfactors(uint64_t perm, uint64_t inverse, int currentsize, int finalsize, int maxpatternsize, hashdb &patterncomplements, hashdb &prefixmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, int prevcount, int* prevextensions) {
  int actualcounts[currentsize + 1];
  int newextensions[currentsize + 1];
  uint64_t newinverse = setdigit(inverse, currentsize, currentsize); // inverse of the extended permutation
  if (currentsize < finalsize) {
    for (int i = currentsize; i >= 0; i--) {
      if (i < currentsize) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * currentsize));
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      int countusingfinaltwo = 0;
      checkpatterns(extendedperm, newinverse, 0, 0, currentsize + 1, maxpatternsize, 0, patterncomplements, prefixmap, countusingfinaltwo);
      //tempcount++;
      //if (tempcount == 1000000L) cout<<tempcount<<endl;;
      int indexignoringsecondtofinal = i;
      if (getdigit(newinverse, currentsize) > getdigit(newinverse, currentsize - 1)) indexignoringsecondtofinal--; 
      int actualcount = countusingfinaltwo + prevcount + prevextensions[indexignoringsecondtofinal];
      increasetally(tally[currentsize + 1], actualcount);
      //completelist[currentsize + 1][permtonum(extendedperm, currentsize + 1)] = actualcount;
      newextensions[i] = countusingfinaltwo + prevextensions[indexignoringsecondtofinal];
      actualcounts[i] = actualcount;
      //displayperm(extendedperm);
      //cout<<"prevcount countusingfinaltwo prevextensionsportion "<<prevcount<<" "<<countusingfinaltwo<<" "<<prevextensions[indexignoringsecondtofinal]<<endl;
    }  
  }
  
  newinverse = setdigit(inverse, currentsize, currentsize); // inverse of the extended permutation
  if (currentsize < finalsize - 1) {
    for (int i = currentsize; i >= 0; i--) {
      if (i < currentsize) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * currentsize));
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations_brute_usingbothaddfactors(extendedperm, newinverse, currentsize + 1, finalsize, maxpatternsize, patterncomplements, prefixmap, tally, completelist, actualcounts[i], newextensions);
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
  timestamp_t start_time = get_timestamp();
  cout<<"max size: "<<maxpermsize<<endl;
  if (!USESECONDADDFACTOR) {
    buildpermutations_brute(0L, 0L, 1, maxpermsize, maxpatternsize, maxpermsize, patterncomplements, prefixmap, tally, completelist, 0);
  } else {
    int temparray[2] = {0, 0};
    buildpermutations_brute_usingbothaddfactors(0L, 0L, 1, maxpermsize, maxpatternsize, patterncomplements, prefixmap, tally, completelist, 0, temparray);
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
// Updates tally and completelist using P_0(perm)
static void Pcount(uint64_t perm, int length, int maxpatternsize, int maxpermsize, const hashdb &patternset, const hashdb &prefixmap, hashmap &Phashmap, uint64_t prevP0, uint64_t prevP1,  vector < vector < int > > &tally, vector < vector < int > > &completelist) {
  stat1++;
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
	stat2++;
      }
    }
    oldPvalpos++;
    if (USEPREFIXMAP) {
      extendnormalizetop(perm, inverse, length, length - 1 - i, prefixentry, seenpos);
      if (!prefixmap.contains(prefixentry)) {
	//displayperm(perm);
	//cout<<oldPvalpos<<endl;
	for (int j = oldPvalpos; j <= Pvalpos; j++) {
	  oldPvals[j] = 0;  // TECHNICALLY THIS MAKES US NOT O(N!) TIME LIKE WE SHOULD BE. BUT NOT WORTH CHANGING IN PRACTICE. BECAUSE THE CACHE MISSES FROM GETTING ACTUAL PVALS IS WHAT DRIVES THE PERFORMANCE, NOT SOME LOOP FILLING THINGS IN WITH ZEROS. STILL, CHANGE THIS IF I GET A CHANCE
	}
	// cout<<"broke"<<endl;
	break;
      } else {
	//if (i < length - 1)  cout<<"Here"<<endl;
      }
    }
  }
  //for (int i = 0; i <= Pvalpos; i++) cout<<oldPvals[i]<<" ";
  //cout<<endl;

  while (Pvalpos >= 0) {
    Pvals[Pvalpos] = Pvals[Pvalpos + 1] + oldPvals[Pvalpos]; // use recurrence to get actual Pvals
    Pvalpos--;
  }
  if (length != maxpermsize) {
    setPvals(perm, Pvals, Phashmap);
  }
    // CAN SAVE CONSIDERABLE TIME (FOR ONE-PATTERN SET) IF I MAKE NEXT TWO LINES OPTIONAL
  //int index = permtonum(perm, length);
  //completelist[length][index] = Pvals[0];
  increasetally(tally[length], Pvals[0]);
  //  displayperm(perm);
  //cout<<"num hits: "<<Pvals[0]<<endl;
}


// constructs the permutations of size finalsize. Then passes on paramaters to Pcount, in order to find number of patterns appearing in each permutation
// Note first call to this function should be with currentsize = 1
void buildpermutations(uint64_t perm, int currentsize, int finalsize, int maxpatternsize, int maxpermsize, hashdb &patternset, hashdb &prefixmap, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, uint64_t prevP0, uint64_t *cachedP1s) {
  if (USEOLDP1 && currentsize == finalsize - 2) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      cachedP1s[i] = getPval(extendedperm, 1, Phashmap);
    }
  }
  uint64_t currentP0 = 0;
  if (USEOLDP0 && currentsize > 1) currentP0 = getPval(perm, 0, Phashmap);
  bool nseen = false;
  for (int i = 0; i < currentsize + 1; i++) {
    if (i > 0 && getdigit(perm, i - 1) == currentsize - 1) nseen = true;
    uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
    if (currentsize < finalsize - 1) {
      buildpermutations(extendedperm, currentsize + 1, finalsize, maxpatternsize,  maxpermsize, patternset, prefixmap, Phashmap, tally, completelist, currentP0, cachedP1s);
    } else {
      uint64_t prevP1 = cachedP1s[i];
      if (nseen) prevP1 = cachedP1s[i - 1];
      Pcount(extendedperm, currentsize + 1, maxpatternsize, maxpermsize, patternset, prefixmap, Phashmap, currentP0, prevP1, tally, completelist);
      //tempcount++;
      //if (tempcount == 1000000L) cout<<tempcount<<endl;;
    }
  }
}

// Fills in tally, complete list, and Pvals for all permutation in S_{<= finalsize}
// Note: patterns in patternset required to be in S_{>1}
void createPmap(uint64_t finalsize, hashdb &patternset, int maxpatternsize, timestamp_t start_time, hashmap &Phashmap, vector < vector < int > > &tally, vector < vector < int > > &completelist, bool verbose) {
  unsigned short temp[maxpatternsize + 1];
  for (int i = 0; i < maxpatternsize + 1; i++) temp[i] = 0;
  setPvals(0L, temp, Phashmap); // fill in Pvals to be 0 for

  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap);

  for (int i = 2; i <= finalsize; i++) {
    uint64_t cachedP1s[finalsize];
    for (int x = 0; x < finalsize; x++) cachedP1s[x] = 0; // initialize to zero used in first call to buildpermutations when i = 2
    buildpermutations(0L, 1, i, maxpatternsize, finalsize, patternset, prefixmap, Phashmap, tally, completelist, 0L, cachedP1s);
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
  stat1 = 0;
  stat2 = 0;
  if (USEBRUTE)  start_brute(maxpatternsize, maxpermsize, patternset, tally, completelist);
  else createPmap(maxpermsize, patternset, maxpatternsize, current_time, Phashmap, tally, completelist, verbose);  //cout<<"Average number of lookups per permutation checked: "<<(double)stat2 / (double)stat1<<endl;
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
void countpatterns(string patternlist, int maxpermsize, vector < vector <int> > & tally, vector < vector < int > > &completelist, bool verbose) {
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
  stat1 = 0;
  stat2 = 0;
  if (USEBRUTE)  start_brute(maxpatternsize, maxpermsize, patternset, tally, completelist);
  else createPmap(maxpermsize, patternset, maxpatternsize, current_time, Phashmap, tally, completelist, verbose);
  cout<<"Average number of lookups per permutation checked: "<<(double)stat2 / (double)stat1<<endl;
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
