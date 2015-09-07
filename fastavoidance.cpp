#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>      /* printf, scanf, puts, NULL */
#include <cstdlib>
#include <cstring>
#include <ctime> 
#include <iostream>
#include <fstream>
#include <queue>
#include <sys/time.h>
#include <vector>
#include "hashdb.h"
#include "fastavoidance.h"
#include "perm.h"
using namespace std;

// Conventions for comments: We use all the conventions from our paper in our comments. Thus the comments disagree with the code in that they do not zero-index the values or positions of a permutation

#define USEBITHACK 1 // 1 to use a bithack inspired by permlab; implemented both for brute force and non-brute force algs.
#define USEPREFIXMAP 1 // only has meaning in non-brute force algorithm. 1 to use the trick which checks for each i-prefix of w whether it is order-isomorphic to an i-prefix of some \pi \in \Pi
#define SINGLEPATTERNOPT 0 // 1 if you want brute-force algorithm to test for each pattern separately rather than use hash table of pattern prefixes to check for all patterns at once whether a subsequence is order isomorphic to any pattern prefixes. Gets some speedup for single-pattern case.
#define USEBRUTE 1 // whether to use brute-force algorithm
#define VERBOSE 0 // whether to be verbose. normally should be false
#define GETSTAT 1 // whether or not to collect statistics -- slows things down a bit. Only used in function run_interior_experiment, and is optionally used for countavoidersfromfile when verbose argument is given to function. Is not implemented for SINGLEPATTERNOPT

// IN CASE OF GETSTAT:
// For brute force variants:
// stat1 counts over all permutations in S_n how many subsequences of length at least three we look at in the course of the algorithm -- unless we are using USEBITHACK, in which case stat1 = stat4
// stat2 is the number of times stat1 is incremented for permutations that end up being avoiders
// stat3 counts total number of avoiders in S_n
// stat4 counts how many times checkpatterns is called in total
// For non-brute-force variants:
// stat1 is total calls to isavoider
// stat2 is number of i-prefixes we look at in total in isavoider calls. if USEBITHACK, excludes 1- and 2-prefixes.
// stat3 counts total number of avoiders in S_n
// stat4 counts number of times stat2 is incremented for permutations that end up being avoiders.
static unsigned long long stat1 = 0, stat2 = 0, stat3 = 0, stat4 = 0; // only used for testing purposes in getstat
static bool countstat1; // is set to true for permutations of size n.


// Get pos-th bit from bitmap. Indexing starts at 0
inline uint64_t getbit(uint64_t word, int pos) {
  return (word >> pos) & 1;
}

// Set pos-th bit from bitmap. Indexing starts at 0
inline uint64_t setbit(uint64_t word, uint64_t pos, uint64_t val) { // val is 1 or 0
  return (word & (~(0L) - (1<<pos))) + (1<<pos)*val;
}

// Shift the bits in positions pos, pos+1, ... to the right and insert bit with value val in position pos. Indexing starts at 0
// Note that due to bit-shift setup, is not capable to inserting in final position ==> Prereq: pos < 63
inline uint64_t insertbit(uint64_t word, uint64_t pos, uint64_t val) { // val is 1 or 0
  return (word & ((1 << pos) - 1) ) + ((word >> pos) << (pos + 1)) + val * (1 << pos);
}


// Note: to prevent a function from inlining for debugging, use this before the function decleration:
//__attribute__((noinline))


// Recursively checks whether perm contains a pattern from patternset. However, prefixmap and patternset
// Requires that all patterns in patternset are length currentpatternlength + numlettersleft
// currentpatterncomplement stores complement of normalization of permutation subsequence already examined
// largestletterused tells us the value of the largest letter icnluded so far in the permutation-subsequence being examined
// seenpos is a bitmap used to efficiently update currentpatterncomplement as we make the subsequence-being looked at one longer
// prefixmap contains the complements of each prefix of each pattern in \Pi.
// however, if SINGLEPATTERNOPT, then prefixes is used in place of prefixmap; prefixes[i] gives us the complement of the i-th normalized prefix of the pattern we are currently searching for.
// Note that prefixmap contains complements of normalized prefixes of patterns 
// Returns true if permutation subsequence cannot be completed to get a pattern, false otherwise
static bool checkpatterns(uint64_t perm, uint64_t inverse, uint64_t currentpatterncomplement, int currentpatternlength, int largestletterused, int numlettersleft, uint32_t seenpos, const hashdb &prefixmap, const vector < uint64_t > & prefixes) {
  if (GETSTAT) stat4++;
  if (GETSTAT && (currentpatternlength >= 3 || !USEBITHACK) && countstat1) stat1++;
  if (!SINGLEPATTERNOPT && currentpatternlength > 1 && !prefixmap.contains(currentpatterncomplement)) return true; 
  if (numlettersleft == 0) return false; // At this point, we have found a pattern
  for (int i = largestletterused - 1; i >= numlettersleft - 1; i--) { // looking at candidates to add to current permutation subsequence
    if (currentpatternlength == 0 && i < largestletterused - 1) return true; // because of how we build candidates for S_n(pi), we can stop here; every permutation this function operates on will have only patterns using n.
    if (USEBITHACK && currentpatternlength == 1 && i < largestletterused - 1) return true; // bithack tells us only need to worry about patterns using both of n-1 and n

    // Similarly to as in extendnormalizetop (defined in perm.cpp), we will build the complement of the normalization of the new permutation-subsequence (with the new letter added to it)
    int oldpos = getdigit(inverse, i);
    int newpos = 0;
    if (oldpos != 0){  
      uint32_t temp = seenpos << (32 - oldpos); // Note: shifting by 32 is ill-defined, which is why we explicitly eliminate digit = 0 case.
      newpos = __builtin_popcount(temp);
    }
    uint64_t newpattern = setdigit(addpos(currentpatterncomplement, newpos), newpos, currentpatternlength);

    // Continue down recursion of either (1) we are not using SINGLEPATTERNOPT, or (2) our newpattern is a valid complement of some prefix of the pattern we are interested in
    if (!SINGLEPATTERNOPT || currentpatternlength == 0 || prefixes[currentpatternlength + 1] == newpattern) { 
      if (checkpatterns(perm, inverse, newpattern, currentpatternlength + 1, i, numlettersleft - 1, seenpos | (1 << oldpos), prefixmap, prefixes) == false) return false;
      // update seenpos and currentpatternelength to be correct as arguments for next step in recursion
    }
  }
  return true;
}


// Checks using brute force algorithm whether permutation is avoider
static bool isavoider_brute(uint64_t perm, uint64_t inverse, int maxavoidsize, int length, const hashdb &prefixmap, const vector < vector < uint64_t > > & prefixes) {
  if (SINGLEPATTERNOPT) {
    for (int i = 0; i < prefixes.size(); i++) {
      if (!checkpatterns(perm, inverse, 0, 0, length, maxavoidsize, 0, prefixmap, prefixes[i])) return false;
    }
    return true;
  } else {
    uint64_t oldstat1 = stat1;
    bool answer =  checkpatterns(perm, inverse, 0, 0, length, maxavoidsize, 0, prefixmap, prefixes[0]);
    if (GETSTAT && answer == true && countstat1) stat2 = stat2 + stat1 - oldstat1;
    return answer;
  }
}

void buildavoiders_brute_helper(uint64_t perm, uint64_t inverse, uint64_t length, uint32_t bitmap, const hashdb &prefixmap, const vector < vector < uint64_t > > & prefixes,int maxavoidsize, int maxsize,  vector < vector < uint64_t > > &avoidervector, vector < uint64_t > &numavoiders, bool justcount) {
  uint64_t newinverse = setdigit(inverse, length, length); // inverse of the extended permutation
  uint64_t newinverses[length+1]; // inverses of each of the extended perms
  uint64_t newperms[length+1]; // each of the extended perms
  for (int i = length; i >= 0; i--) {
    // need to increment newinverse[perm[i]], decrement newinverse[length]
    if (i < length) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * length));
    newinverses[i] = newinverse;
    uint64_t extendedperm = setdigit(addpos(perm, i), i, length); // insert length in i-th position (remember, values are indexed starting at 0)
    newperms[i] = extendedperm;
    if (!USEBITHACK || getbit(bitmap, i) == 1) { // If we are using bithack, then we only bother extending perm by inserting value length in i-th position if the bitmap tells tells us the result is a potential avoider
      if (GETSTAT) countstat1 = false;
      if (GETSTAT && maxsize == length + 1) countstat1 = true;
      if (isavoider_brute(extendedperm, newinverse, maxavoidsize, length + 1, prefixmap, prefixes)) { // if extended permutation is avoider
	if (!justcount) avoidervector[length + 1].push_back(extendedperm);
	else numavoiders[length + 1]++;
      } else {
	newperms[i] = -1; // signifies that we should NOT go further down recursion, for when not using bithack
	if (USEBITHACK) bitmap = setbit(bitmap, i, 0); // keep track of which insertion positions resulted in an avoider 
      }
    }
  }
  uint32_t newmap = bitmap;
  if (length + 1 < maxsize) {
    for (int i = length; i >= 0; i--) {
      if (!USEBITHACK || getbit(bitmap, i) == 1) {
	if (USEBITHACK) newmap = insertbit(bitmap, i + 1, 1);
	if (USEBITHACK || newperms[i] != -1) buildavoiders_brute_helper(newperms[i], newinverses[i], length + 1, newmap, prefixmap, prefixes, maxavoidsize, maxsize,  avoidervector, numavoiders, justcount);
	//bitmaps.push(insertbit(bitmap, i + 1, 1)); // using which insertion positions resulted in an avoider, build bitmap for each new avoider
      }
    }
  }
}

void buildavoiders_brute(const hashdb &patternset, int maxavoidsize, int maxsize,  vector < vector < uint64_t > > &avoidervector, vector < uint64_t > &numavoiders, bool justcount, uint64_t plannedavoidsetsize) {
  stat1 = 0;
  stat2 = 0;
  stat3 = 0;
  stat4 = 0;
  if (VERBOSE && USEBITHACK) cout<<"Using Brute Force Algorithm, Based on PermLab's"<<endl;
  if (VERBOSE && !USEBITHACK) cout<<"Using naive brute force algorithm"<<endl;
  if (!justcount) avoidervector.resize(maxsize + 1);
  else numavoiders.resize(maxsize + 1);
  
  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap); // defined in perm.cpp
  vector < vector < uint64_t > > prefixes;
  prefixes.resize(patternset.getsize());

  vector <unsigned long long> patterns;
  //hashdb patterncomplements(1<<3);
  patternset.getvals(patterns);
  for (int i = 0; i < patterns.size(); i++) {
    uint64_t perm = patterns[i];
    int length = getmaxdigit(perm) + 1;
    //patterncomplements.add(getcomplement(perm, length));

    prefixes[i].resize(length + 1);
    uint64_t entry = 0;
    uint64_t inverse = getinverse(perm, length);
    uint32_t seenpos = 0; // bit map of which positions we've seen so far
    for (int j = 0; j < length; j++) {
      extendnormalizetop(perm, inverse, length, j, entry, seenpos); // defined in perm.cpp
      prefixes[i][j + 1] = entry;
    }
  }

  buildavoiders_brute_helper(0L, 0L, 0, 1, prefixmap, prefixes, maxavoidsize, maxsize, avoidervector, numavoiders, justcount);
}


// Detects if perm is in S_{length}(patternset), where maxavoidsize is the length of the longest pattern in patternset.
// Prerequisite:
// -- The subset of S_{n-1} contained in avoidset is exactly S_{n-1}(patternset)
// -- If USEBITHACK, then we also require that the user has already verified that all patterns in perm use both n and n-1.
// -- If USEPREFIXMAP, prefixmap is used to determine whether prefixes of perm are order isomorphic to prefixes of permutations in patternset
// prefixmap actually contains the complements of the normalizations of the prefixes of the patterns.
static bool isavoider(uint64_t perm, uint64_t inverse, int maxavoidsize, int length, const hashdb &avoidset, const hashdb &patternset, const hashdb &prefixmap) {
  if (GETSTAT) stat1++;
  if (length <= maxavoidsize && patternset.contains(perm)) { // if perm is an offending pattern
      return false;
  }
  uint32_t seenpos = 0; // will be a bitmap used for prefix creation
  uint64_t prefixentry = 0; // will contains the coplement of the normalization of the prefix of perm we're currently looking at
  uint64_t currentperm = perm;
  if (length > 1) { // don't deal with permutations of size zero
    for (int i = length - 1; i >= 0 && i >= length - maxavoidsize - 1; i--) { // for i ranging from the largest-valued letter in perm to the (maxavoidsize + 1)-th largest-valued letter in perm
      // Note: the length - 1 case is only for setting up the prefix stuff correctly, and because we compute perm \downarrow_2 from perm \downarrow_1
      if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
        currentperm = addpos(currentperm, getdigit(inverse, i + 1));
	currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
      }
      currentperm = killpos(currentperm, getdigit(inverse, i)); // now currentperm is perm, except with the letter of value i removed, and the permutation normalized to be on the letters 0,...,(length - 2)
      if (GETSTAT && (!USEBITHACK || i < length - 2)) stat2++;
      if ((!USEBITHACK || i < length - 2) && !avoidset.contains(currentperm)) { // Check if the prefix is in S_{length - 1}(patternset)
	return false; // found a subword not avoiding the patterns
      }

      if (USEPREFIXMAP) {
	extendnormalizetop(perm, inverse, length, length - 1 - i, prefixentry, seenpos); // defined in perm.cpp
	if (i < length - 1 && !prefixmap.contains(prefixentry)) return true; 
	if (i == length - maxavoidsize) return false; // in this case, the prefix must actually be a pattern
      }
    }
  }
  return true;
}


// Builds the permutations in $S_1, ..., S_maxsize$ avoiding the
// patterns in patternset, the longest of which is length
// maxavoidsize. If justcount, does nothing with avoidervector but
// makes nuavoiders[i] be the number of avoiders in S_i (for i >
// 0). If !justcount, does nothing with numavoiders, but makes
// avoidervector contain a vector of all permutations in S_i in
// avoidervector[i] (for i > 0). plannedavoisetsize should be large if
// we expect a major computation and small otherwise (dictates how
// much memory is initially allocated to data structures at start of
// algorithm.).
// Note: Patternset patterns required to be size >= 2
void buildavoiders(const hashdb &patternset, int maxavoidsize, int maxsize,  vector < vector < uint64_t > > &avoidervector, vector < uint64_t > &numavoiders, bool justcount, uint64_t plannedavoidsetsize) {
  if (!justcount) avoidervector.resize(maxsize + 1);
  else numavoiders.resize(maxsize + 1);

  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap);
  if (VERBOSE) cout<<"Planned size: "<<plannedavoidsetsize<<endl;
  hashdb avoidset = hashdb(plannedavoidsetsize); // hash table containing avoiders of all sizes
  uint64_t startperm = 0;
  avoidset.add(startperm); // identity in S_1
  if (!justcount) avoidervector[1].push_back(startperm);
  else numavoiders[1] = 1;
  
  std::queue<unsigned long long> avoiderstoextend; // queue of avoiders built so far.
  // when we find an avoider, we will add it to this queue. We will
  // then later take it out of the queue and use it to generate
  // options for avoiders of length one larger
  std::queue<unsigned long long> bitmaps; // Contains a bitmap associated with each avoider in avoiderstoextendd (if USEBITHACK)
  // The bitmap for a permutation w \in S_n has a 1 in position i iff
  // inserting letter (n + 1) in position i of w would result in
  // permutation w' such that if you removed the letter n from w', the
  // result would be a patternset-avoiding word. Thus when checking
  // whether w' is an avoider, we do not have to explicitly check for
  // this property, preventing a potential cache-miss.
  
  avoiderstoextend.push(startperm);
  if (USEBITHACK) bitmaps.push(3L);

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
    uint64_t bitmap = bitmaps.front();
    avoiderstoextend.pop();
    if (USEBITHACK) bitmaps.pop();
    numleftcurrentlength--;
    uint64_t inverse = getinverse(perm, currentlength);
    uint64_t newinverse = setdigit(inverse, currentlength, currentlength); // inverse of the extended permutation
    for (int i = currentlength; i >= 0; i--) {
      // need to increment newinverse[perm[i]], decrement newinverse[currentlength]
      if (i < currentlength) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * currentlength));
      if (!USEBITHACK || getbit(bitmap, i) == 1) { // If we are using bithack, then we only bother extending perm by inserting value currentlength in i-th position if the bitmap tells tells us the result is a potential avoider
	uint64_t extendedperm = setdigit(addpos(perm, i), i, currentlength); // insert currentlength in i-th position (remember, values are indexed starting at 0, so results in permutation in S_currentlength)
	uint64_t tempstat2 = stat2;
	if (isavoider(extendedperm, newinverse, maxavoidsize, currentlength + 1, avoidset, patternset, prefixmap)) { // if extended permutation is avoider
	  if (GETSTAT) stat4 += stat2 - tempstat2;
	  if (!justcount) avoidervector[currentlength + 1].push_back(extendedperm);
	  else numavoiders[currentlength + 1]++;
	  if (currentlength + 1 < maxsize) {
	    avoiderstoextend.push(extendedperm);
	    avoidset.add(extendedperm);
	    numnextlength++; 
	  }
	} else {
	  if (USEBITHACK) bitmap = setbit(bitmap, i, 0); // keep track of which insertion positions resulted in an avoider
	}
      }
    }
    if (USEBITHACK && currentlength + 1 < maxsize) {
      for (int i = currentlength; i >= 0; i--) {
	if (getbit(bitmap, i) == 1) {
	  bitmaps.push(insertbit(bitmap, i + 1, 1)); // using which insertion positions resulted in an avoider, build bitmap for each new avoider
	}
      }
    }
  }
}



// Example:
// string patternlist = "1234 3214"; // space separated list of patterns; need not be same sizes; must be in S_{<10}
// vector < vector < uint64_t > > avoidervector;
// buildavoidersfrompatternlist(patternlist, 10, avoidervector); // now avoidervector contains S_n(patternlist) stored in avoidervector[n] for 0 < n < 11
void buildavoidersfrompatternlist(string patternlist, int maxpermsize, vector < vector < uint64_t > > &avoidervector) {
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);
  vector < uint64_t > numavoiders;
  if (USEBRUTE) buildavoiders_brute(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, false, (1L << 10));
  else buildavoiders(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, false, (1L << 10)); // for large cases, make last argument much larger
}

// just for use when GETSTATIS TRUE
uint64_t getstat1() {
  return stat1;
}
uint64_t getstat2() {
  return stat2;
}
uint64_t getstat3() {
  return stat3;
}
uint64_t getstat4() {
  return stat4;
}

// Useful for analyzing how much time is spent where. Should be run with GETSTAT on
double run_interior_experiment(string patternlist, int maxpermsize) {
  stat1 = 0;
  stat2 = 0;
  stat3 = 0;
  stat4 = 0;
  timestamp_t start_time = get_timestamp();
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);
  vector < uint64_t > numavoiders;
  vector < vector < uint64_t > > avoidervector;
  if (VERBOSE) cout<<"Effective pattern size "<<maxpatternsize<<endl;
  if (USEBRUTE) buildavoiders_brute(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true, (1L << 10));
  else buildavoiders(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true, (1L << 10)); // for large cases, could make last argument much larger. But we will not bother in our tests in the paper.
  stat3 = numavoiders[maxpermsize];
  timestamp_t end_time = get_timestamp();
  return (end_time - start_time)/1000000.0L;
}

// Example:
// string patternlist = "1234 3214"; // space separated list of patterns; need not be same sizes; must be in S_{<10}
// vector < uint64_t > numavoiders;
// buildavoidersfrompatternlist(patternlist, 10, numavoiders); // now avoidervector contains |S_n(patternlist)| stored in numavoiders[n] for 0 < n < 11.
void countavoidersfrompatternlist(string patternlist, int maxpermsize, vector < uint64_t > &numavoiders) {
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);
  vector < vector < uint64_t > > avoidervector;
  if (VERBOSE) cout<<"Effective pattern size "<<maxpatternsize<<endl;
  if (USEBRUTE) buildavoiders_brute(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true, (1L << 10));
  else buildavoiders(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true, (1L << 10)); // for large cases, make last argument much larger
}

// Inputs file stream containing string list of patterns on each line. 
// e.g., infile might contain
// 1234 312
// 231 4132 312
// Outputs file stream alternating every other line
// (1) #<list of patterns>
// (2) |S_1(list of patterns)| |S_2(list of patterns)| ... |S_maxpermsize(list of patterns)|
// e.g., output file might contain
// #1234 312
// 0 2 5 13 31 66 127 225 373 586 
// #231 4132 312
// 0 2 4 8 16 32 64 128 256 512
void countavoidersfromfile(ifstream &infile, ofstream &outfile, int maxpermsize, bool verbose) {
  string line;
  while (getline(infile, line)) {
    outfile<<"#"<<line<<endl;
    if (verbose) cout<<line<<endl;
    vector < uint64_t > numavoiders;
    timestamp_t start_time = get_timestamp();
    countavoidersfrompatternlist(line, maxpermsize, numavoiders);
    timestamp_t end_time = get_timestamp();
    for (int i = 1; i < numavoiders.size(); i++) {
      if (verbose) cout<<numavoiders[i]<<" ";
      outfile<<numavoiders[i]<<" ";
    }
    if (verbose) cout<<endl;
    outfile<<endl;
    if (verbose) cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
    stat3 = numavoiders[maxpermsize];
    if (verbose && GETSTAT && !USEBRUTE) cout<<(double)stat2/(double)stat1<<" prefixes looked at on average per call to isavoider"<<endl;
    if (verbose && GETSTAT && !USEBRUTE) cout<<(double)stat2/(double)stat3<<" prefixes looked at on average per actual avoider"<<endl;
  }
  return;
}
