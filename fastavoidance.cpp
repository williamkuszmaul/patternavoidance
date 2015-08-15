// simplest make is g++ -O2 -std=c++11  fastavoidance.cpp  -g -o fastavoidance
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

// 1 to use a bithack inspired by permlab. Upshot: Code runs about a third faster. Downshot: Code is around 1/4 less memory efficient
#define USEBITHACK 1
#define USEPREFIXMAP 1

// Input: perm, inverse (which needs to be right in pos length - index - 1), perm length, index, complement of normalization of index - 1 largest letters in perm, a bitmap which should start off at zero for index = 0
// Output: bitmap is updated. answer is updated to be complement of normalization of index largest letters in perm
void extendnormalizetop(uint64_t perm, uint64_t inverse, int length, int index, uint64_t &answer, uint32_t & seenpos) {
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

void addprefixeshelper(uint64_t perm, int length, hashdb &table) {
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
void addprefixes(const hashdb &permset, hashdb &table) {
  vector <unsigned long long> patterns;
  permset.getvals(patterns);
  for (int i = 0; i < patterns.size(); i++) {
    uint64_t perm = patterns[i];
    int length = getmaxdigit(perm) + 1;
    displayperm(perm,length);
    addprefixeshelper(perm, length, table);
  }
}

unsigned long long stat1 = 0, stat2 = 0, stat3 = 0;

// THIS IS HOW YOU KEEP A FUNCTION FROM INLINING!
__attribute__((noinline)) uint64_t getinversecover(uint64_t perm, size_t length) {
  //__attribute__((noinline));
  return getinverse(perm, length);
}

// Detects if perm is in S_{length}(patternset), where maxavoidsize is the length of the longest pattern in patternset.
// Prerequisite:
// -- The subset of S_{n-1} contained in avoidset is exactly S_{n-1}(patternset)
// -- If USEBITHACK, then we also require that the user has already verified that all patterns in perm use both n and n-1.
static bool isavoider(uint64_t perm, uint64_t inverse, int maxavoidsize, int length, const hashdb &avoidset, const hashdb &patternset, const hashdb &prefixmap) {
  stat1++;
  if (length <= maxavoidsize && patternset.contains(perm)) { // if perm is an offending pattern
      return false;
  }
  uint32_t seenpos = 0;
  uint64_t prefixentry = 0;
  uint64_t currentperm = perm;
  if (length > 1) { // don't deal with permutations of size zero
    for (int i = length - 1; i >= 0 && i >= length - maxavoidsize - 1; i--) { // for i ranging from the largest-valued letter in perm to the (maxavoidsize + 1)-th largest-valued letter in perm
      if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
        currentperm = addpos(currentperm, getdigit(inverse, i + 1));
	currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
      }
      currentperm = killpos(currentperm, getdigit(inverse, i)); // now currentperm is perm, except with the letter of value i removed, and the permutation normalized to be on the letters 0,...,(length - 2)
      //displayperm(currentperm)
      if (!USEBITHACK || i < length - 2) stat2++;
      // NEVER ACTUALLY HAVE TO DO I = LENGTH - 1
      if ((!USEBITHACK || i < length - 2) && !avoidset.contains(currentperm)) { // Check if this sub-permutation of perm is in S_{length - 1}(patternset)
	return false; // found a subword not avoiding the patterns
      }

      if (USEPREFIXMAP) {
	extendnormalizetop(perm, inverse, length, length - 1 - i, prefixentry, seenpos);
	if (i < length - 1 && !prefixmap.contains(prefixentry)) return true;
	if (i == length - maxavoidsize) return false;
      }
      // The next two lines are special for when there is a single pattern and it is the identity. They reduce cache misses significantly. For |pi|=5, this gets roughly a x2 speedup
      // Can this hack be extended efficiently to arbitrary choices of pattern-sets?
      // if (i < length - 1 && getdigit(inverse, i) > getdigit(inverse, i + 1)) return true;
      // if (i == length - maxavoidsize) return false;
    }
  }
  return true;
}

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
  // this property, preventing a cache-miss.
  
  avoiderstoextend.push(startperm);
  if (USEBITHACK) bitmaps.push(3L);

  int currentlength = 1; // maintain as length of next permutation to be popped from avoiderstoextend
  int numleftcurrentlength = 1; // number of permutations left in avoiderstoextend until we have to increment currentlength
  int numnextlength = 0; // number of permutations of size currentlength + 1 in avoiderstoextend

  while (avoiderstoextend.size() > 0) {
    if (numleftcurrentlength == 0) {
      //cout<<"Finished n = "<<currentlength + 1<<" with "<<numavoiders[currentlength + 1]<<" avoiders"<<endl; // MEMORY EFFICIENCY NOTE: I TEMPORARILY LOOSE A FACTOR OF 1/3 WHEN I RESIZE
      numleftcurrentlength = numnextlength;
      numnextlength = 0;
      currentlength++;
    }
    uint64_t perm = avoiderstoextend.front();
    uint64_t bitmap = bitmaps.front();
    avoiderstoextend.pop();
    if (USEBITHACK) bitmaps.pop();
    numleftcurrentlength--;
    uint64_t inverse = getinversecover(perm, currentlength);
    uint64_t newinverse = setdigit(inverse, currentlength, currentlength); // inverse of the extended permutation
    for (int i = currentlength; i >= 0; i--) {
      // need to increment newinverse[perm[i]], decrement newinverse[currentlength]
      stat3++;
      if (i < currentlength) newinverse = newinverse + (1L << (4 * getdigit(perm, i))) - (1L << (4 * currentlength));
      if (!USEBITHACK || getbit(bitmap, i) == 1) { // If we are using bithack, then we only bother extending perm by inserting value currentlength in i-th position if the bitmap tells tells us the result is a potential avoider
	uint64_t extendedperm = setdigit(addpos(perm, i), i, currentlength); // insert currentlength in i-th position (remember, values are indexed starting at 0)
	if (isavoider(extendedperm, newinverse, maxavoidsize, currentlength + 1, avoidset, patternset, prefixmap)) { // if extended permutation is avoider
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
  buildavoiders(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, false, (1L << 10)); // for large cases, make last argument much larger!
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
  buildavoiders(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true, (1L << 10)); // for large cases, make last argument much larger!
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
  }
  cout<<(double)stat2/(double)stat3<<endl;
  return;
}
