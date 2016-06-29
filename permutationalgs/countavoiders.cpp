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

#ifndef _NO_CILK
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_list.h>
#endif

#include <queue>
#include <sys/time.h>
#include <vector>
#include "countavoiders.h"
#include "hashdb.h"
#include "hashmap.h"
#include "perm.h"
#include "permutilities.h"
using namespace std;

// Conventions for comments: We use all the conventions from our paper in our comments. Thus the comments disagree with the code in that they do not zero-ind ex the values or positions of a permutation

#define USEBITHACK 1
// 1 to use a bithack inspired by permlab; implemented ONLY for brute force. (is default for non-brute force algs).
#define USEPREFIXMAP 1 // only has meaning in non-brute force algorithm. 1 to use the trick which checks for each i-prefix of w whether it is order-isomorphic to an i-prefix of some \pi \in \Pi
#define SINGLEPATTERNOPT 1 // 1 if you want brute-force algorithm to test for each pattern separately rather than use hash table of pattern prefixes to check for all patterns at once whether a subsequence is order isomorphic to any pattern prefixes. Gets some speedup for single-pattern case.
// IF SINGLEPATTERNOPT IS FALSE, THEN FOR BRUTE FORCE ALGORITHM PATTERNS MUST ALL BE SAME SIZES
#define USEBRUTE 0 // whether to use brute-force algorithm
#define VERBOSE 0 // whether to be verbose. normally should be false
#define GETSTAT 0 // whether or not to collect statistics -- slows things down a bit. Only used in function run_interior_experiment, and is optionally used for countavoidersfromfile when verbose argument is given to function. Is currently only implemented for SINGLEPATTERNOPT = false and USEBRUTE = true.
// stat1 counts over all permutations in S_n how many subsequences of length at least two we look at in the course of the algorithm
// stat2 is the number of times stat1 is incremented for permutations that end up being avoiders
// stat3 counts total number of avoiders in S_n
// stat4 counts how many subsequences of length at least two we look at over the course of the algorithm for ALL permutations
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

// returns position of first 1, with indexing starting at 0
inline uint64_t firstonepos(uint64_t word) {
  assert(word != 0);
  return __builtin_ctzl(word);
}


// returns position of last 1, with indexing starting at 0
inline uint64_t lastonepos(uint64_t word) {
  assert(word != 0);
  return 63 - __builtin_clzl(word);
}


// Note: to prevent a function from inlining for debugging, use this before the function decleration:
//__attribute__((noinline))


// Recursively checks whether perm contains a pattern from patternset.
// Requires that all patterns in patternset are length currentpatternlength + numlettersleft
// currentpatterncomplement stores complement of normalization of permutation subsequence already examined
// largestletterused tells us the value of the largest letter icnluded so far in the permutation-subsequence being examined
// seenpos is a bitmap used to efficiently update currentpatterncomplement as we make the subsequence-being looked at one longer
// prefixmap contains the complements of each prefix of each pattern in \Pi.
// however, if SINGLEPATTERNOPT, then prefixes is used in place of prefixmap; prefixes[i] gives us the complement of the i-th normalized prefix of the pattern we are currently searching for.
// Note that prefixmap contains complements of normalized prefixes of patterns
// Returns true if permutation subsequence cannot be completed to get a pattern, false otherwise
static bool checkpatterns(perm_t perm, perm_t inverse, perm_t currentpatterncomplement, int currentpatternlength, int largestletterused, int numlettersleft, uint32_t seenpos, const hashdb &prefixmap, const vector < perm_t > & prefixes) {
  if (GETSTAT && currentpatternlength > 1) stat4++;
  if (GETSTAT && currentpatternlength > 1 && countstat1) stat1++;
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
    perm_t newpattern = setdigit(addpos(currentpatterncomplement, newpos), newpos, currentpatternlength);

    // Continue down recursion of either (1) we are not using SINGLEPATTERNOPT, or (2) our newpattern is a valid complement of some prefix of the pattern we are interested in
    if (!SINGLEPATTERNOPT || currentpatternlength == 0 || prefixes[currentpatternlength + 1] == newpattern) {
      if (checkpatterns(perm, inverse, newpattern, currentpatternlength + 1, i, numlettersleft - 1, seenpos | (1 << oldpos), prefixmap, prefixes) == false) return false;
      // update seenpos and currentpatternelength to be correct as arguments for next step in recursion
    }
  }
  return true;
}


// Checks using brute force algorithm whether permutation is avoider
static bool isavoider_brute(perm_t perm, perm_t inverse, int maxavoidsize, int length, const hashdb &prefixmap, const vector < vector < perm_t > > & prefixes, vector  < uint64_t > patternlengths) {
  if (SINGLEPATTERNOPT) {
    for (int i = 0; i < prefixes.size(); i++) { // check each pattern individually
      if (!checkpatterns(perm, inverse, 0, 0, length, patternlengths[i], 0, prefixmap, prefixes[i])) return false;
    }
    return true;
  } else { // check all the patterns at once // REQUIRES PATTERNS ALL SAME SIZE
    uint64_t oldstat1 = stat1;
    bool answer =  checkpatterns(perm, inverse, 0, 0, length, maxavoidsize, 0, prefixmap, prefixes[0]);
    if (GETSTAT && answer == true && countstat1) stat2 = stat2 + stat1 - oldstat1;
    return answer;
  }
}


void buildavoiders_brute_helper(perm_t perm, perm_t inverse, uint64_t length, uint64_t bitmap, const hashdb &prefixmap, const vector < vector < perm_t > > & prefixes, vector  < uint64_t > & patternlengths, int maxavoidsize, int maxsize,  vector < cilk::reducer< cilk::op_list_append<perm_t> > > &avoidervector, cilk::reducer< cilk::op_add<uint64_t> > *numavoiders, bool justcount) {

  // We take perm and extend it by adding a new largest letter in each possible position
  perm_t newinverse = setdigit(inverse, length, length); // inverse of the extended permutation
  perm_t newinverses[length+1]; // inverses of each of the extended perms
  perm_t newperms[length+1]; // each of the extended perms
  for (int i = length; i >= 0; i--) {
    // need to increment newinverse[perm[i]], decrement newinverse[length]
    if (i < length)  newinverse = decrementdigit(incrementdigit(newinverse, getdigit(perm, i)), length);
    newinverses[i] = newinverse;
    perm_t extendedperm = setdigit(addpos(perm, i), i, length); // insert length in i-th position (remember, values are indexed starting at 0)
    newperms[i] = extendedperm;
    if (!USEBITHACK || getbit(bitmap, i) == 1) { // If we are using bithack, then we only bother extending perm by inserting value length in i-th position if the bitmap tells tells us the result is a potential avoider
      if (GETSTAT) countstat1 = false;
      if (GETSTAT && maxsize == length + 1) countstat1 = true;
      if (isavoider_brute(extendedperm, newinverse, maxavoidsize, length + 1, prefixmap, prefixes, patternlengths)) { // if extended permutation is avoider
	if (!justcount) (*avoidervector[length + 1]).push_back(extendedperm);
	else *(numavoiders[length + 1]) += 1;
      } else {
	newperms[i] = -1; // signifies that we should NOT go further down recursion, for when not using bithack
	if (USEBITHACK) bitmap = setbit(bitmap, i, 0); // keep track of which insertion positions resulted in an avoider
      }
    }
  }

  // For the extended permutations that succeeded at avoiding the patterns, we recurse on them
  uint64_t newmap = bitmap;
  if (length + 1 < maxsize) {
    for (int i = length; i >= 0; i--) {
      if (!USEBITHACK || getbit(bitmap, i) == 1) {
	if (USEBITHACK) newmap = insertbit(bitmap, i + 1, 1);
	if (USEBITHACK || newperms[i] != -1) {
	  cilk_spawn  buildavoiders_brute_helper(newperms[i], newinverses[i], length + 1, newmap, prefixmap, prefixes, patternlengths, maxavoidsize, maxsize,  avoidervector, numavoiders, justcount);
	}
      }
    }
    cilk_sync;
  }
}

void buildavoiders_brute(const hashdb &patternset, int maxavoidsize, int maxsize,  vector < list < perm_t > > avoidervector, vector < uint64_t > &numavoiders, bool justcount, uint64_t plannedavoidsetsize) {
  stat1 = 0;
  stat2 = 0;
  stat3 = 0;
  stat4 = 0;
  if (VERBOSE && USEBITHACK) cout<<"Using Algorithm Based on PermLab's"<<endl;
  if (VERBOSE && !USEBITHACK) cout<<"Using Naive Brute Force Algorithm"<<endl;
  if (!justcount) avoidervector.resize(maxsize + 1);
  else numavoiders.resize(maxsize + 1);

  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap); // defined in perm.cpp
  vector < vector < perm_t > > prefixes;
  vector  < uint64_t > patternlengths; // the length of each pattern
  prefixes.resize(patternset.getsize());
  patternlengths.resize(patternset.getsize());

  vector <perm_t> patterns;
  //hashdb patterncomplements(1<<3);
  patternset.getvals(patterns);
  for (int i = 0; i < patterns.size(); i++) {
    perm_t perm = patterns[i];
    int length = getmaxdigit(perm) + 1;
    //patterncomplements.add(getcomplement(perm, length));

    prefixes[i].resize(length + 1);
    patternlengths[i] = length;
    perm_t entry = 0;
    perm_t inverse = getinverse(perm, length);
    uint32_t seenpos = 0; // bit map of which positions we've seen so far
    for (int j = 0; j < length; j++) {
      extendnormalizetop(perm, inverse, length, j, entry, seenpos); // defined in perm.cpp
      prefixes[i][j + 1] = entry;
    }
  }

  vector < cilk::reducer< cilk::op_list_append<perm_t> > > avoidervectortemp(maxsize + 1);
  cilk::reducer< cilk::op_add<uint64_t> > *numavoiderstemp = new cilk::reducer< cilk::op_add<uint64_t> >[maxsize + 1];
  buildavoiders_brute_helper(0L, 0L, 0, 1, prefixmap, prefixes, patternlengths, maxavoidsize, maxsize, avoidervectortemp, numavoiderstemp, justcount);
  for (int i = 0; i < maxsize + 1; i++) {
    if (!justcount) avoidervector[i] = avoidervectortemp[i].get_value();
    numavoiders[i] = numavoiderstemp[i].get_value();
  }
  delete[] numavoiderstemp;
}

// Updates inverse to be the inverse of the same permutation, except with length inserted in position pos.
// Only responsible for having accurate values of inverse for final partial letters of new inverse.
// O(partial) time
inline perm_t extendpartialinverse(perm_t inverse, uint64_t pos, uint64_t length, uint64_t partial) {
  perm_t new_inverse = setdigit(inverse, length, pos); // add the new letter to inverse
  for (int i = 1; i < partial && i < length + 1; i++) {
    int letter = length - i;
    int letter_pos = getdigit(inverse, letter);
    if (letter_pos >= pos) {
      new_inverse = setdigit(new_inverse, letter, letter_pos + 1);
    }
  }
  // cout<<"-----------"<<endl;
  // cout<<partial<<endl;
  // cout<<pos<<endl;
  // displayperm(inverse);
  // displayperm(new_inverse);
  return new_inverse;
}


// Implements the depth-first search described in memory-efficient version of dynamic algorithm in paper (see README)
// Moreover, uses bithacks to analyze all the extensions of an avoider at once. (and figure out if they're avoiders)
// The candidates are the permutations which can be extended to obtain potential elements of C^k(perm) \cap S_{\le n}(\Pi)
// candidate_start_pos is position of first candidate in candidates, candidate_end_pos position of last
// avoiderset_read is the hash table C^k(perm\downarrow_1) \cap S_{\le n}(\Pi)
void buildavoiders_dynamic_helper(uint64_t candidate_length, const hashdb &patternset, const hashdb &shiftedprefixmap, int maxavoidsize, int maxsize,  vector < cilk::reducer< cilk::op_list_append<perm_t> > > &avoidervector, cilk::reducer< cilk::op_add<uint64_t> > *numavoiders, bool justcount, const vector <perm_t> &candidates, const vector <perm_t> &candidate_inverses, int64_t candidate_start_pos, int64_t candidate_end_pos, const hashmap &avoidermap_read, const vector <uint64_t> &prevbitmaps) {
  hashmap avoidermap_write(1 << 8, 8);
  vector <perm_t> next_candidates; // the avoiders obtained from current candidates
  vector <perm_t> next_candidate_inverses; // the avoiders obtained from current candidates
  vector <uint64_t> bitmaps;
  
  // Go through the candidates and build for each a bitmap of which positions the candidate can be extended in to get a new avoider
  for (int64_t pos = candidate_start_pos; pos <= candidate_end_pos; pos++) { // signed because end_pos might be -1
    perm_t candidate = candidates[pos];
    perm_t candidate_inverse = candidate_inverses[pos]; // getinverse(candidate, candidate_length);
    uint64_t bitmap = prevbitmaps[pos];
    // based on candidate\downarrow_1, bitmap tells us options for positions we can extend candidate in to
    // get a new avoider.
    uint64_t avoidbits = bitmap; // Will end up saying in which positions candidate can be extended to get a new avoider
    perm_t prefixentry = 0; // will contain the complement of the normalization of the prefix of perm we're currently looking at
    uint32_t seenpos = 0;
    // delete largest letter from candidate to get smaller perm
    perm_t smallerperm = killpos(candidate, getdigit(candidate_inverse, candidate_length - 1));
    // smallerperm will take values of candidate\downarrow_i for i from 2 to maxavoidsize - 1
    // (currently takes the value for i = 1)

    // initialize bitmap stuff needed for USEPREFIXMAP
    if (USEPREFIXMAP) extendnormalizetop(candidate, candidate_inverse, candidate_length, 0, prefixentry, seenpos); // defined in perm.cpp

    for (int killval = candidate_length - 2; killval > (int)candidate_length - 1 - maxavoidsize && killval >= 0; killval--) {
      // add back in digit we deleted a moment ago, but with value one smaller
      smallerperm = addpos(smallerperm, getdigit(candidate_inverse, killval + 1));
      smallerperm = setdigit(smallerperm, getdigit(candidate_inverse, killval + 1), killval);
      // remove the value to be killed
      uint64_t killedpos = getdigit(candidate_inverse, killval);
      smallerperm = killpos(smallerperm, killedpos); // now smallerperm is candidate, except with the letter of value killval removed, and the permutation normalized to be on the letters 0,...,(candidate_length - 2)

      // now patternbits is the bitmap for what positions we can extend smallerperm in to get an avoider
      uint64_t patternbits = *((uint64_t *)(avoidermap_read.getpayload(smallerperm)));
      // We will now translate this into a bitmap telling us possible positions we can extend candidate in to get a new avoider
      // need to copy the killedpos-th bit to killedpos-th position + 1, and slide higher indexed bits each one over to make room
      uint64_t mask1 = ~((1 << (killedpos)) - 1); // all bits after and including killed pos
      uint64_t mask2 = ((1 << (killedpos + 1)) - 1); // all bits before and including killed pos
      patternbits = ((patternbits & mask1) << 1) + (patternbits & mask2);

      // update avoidbits
      avoidbits = avoidbits & patternbits;

      if (USEPREFIXMAP) {
	extendnormalizetop(candidate, candidate_inverse, candidate_length, candidate_length - 1 - killval, prefixentry, seenpos); // defined in perm.cpp
	if (!shiftedprefixmap.contains(prefixentry)) {
	  // Since the largest (candidate_length -  killval) letters of candidate down't form w\downarrow_1 for any pattern w
	  // it follows that any extension of candidate has no pattern using its largest (candidate_length - killval + 1) letters.
	  break;
	}
      }
    }

    // Next we update avoidbits for cases where one of the extensions of candidate are actually a pattern
    if (candidate_length < maxavoidsize) {
      uint64_t tempbits = avoidbits;
      while (tempbits != 0) {
	int i = firstonepos(tempbits); // next position which avoidbits has a one in
	perm_t extendedperm = setdigit(addpos(candidate, i), i, candidate_length); // insert candidate_length in i-th position (remember, values are indexed starting at 0, so results in permutation in S_candidate_length);
	if (patternset.contains(extendedperm)) avoidbits = setbit(avoidbits, i, 0);
	tempbits = tempbits - (1 << i);
      }
    }

    // Now avoidbits correctly tells us which extensions of candidate result in avoiders
    if (candidate_length + 1 < maxsize) avoidermap_write.add(candidate, &avoidbits); // adding even when avoidbits all zeros for now. allows us to not worry about NULLs in lookups

    // If we're at our final size, report the new avoiders found
    if (candidate_length + 1 == maxsize) {
      if (justcount) {
	*(numavoiders[candidate_length + 1]) += __builtin_popcount((uint32_t)avoidbits);
      } else {
	uint64_t tempbits = avoidbits;
	while (tempbits != 0) {
	  int i = lastonepos(tempbits); // next position which avoidbits has a one in
	  tempbits = tempbits - (1 << i);
	  perm_t extendedperm = setdigit(addpos(candidate, i), i, candidate_length); // insert candidate_length in i-th position (remember, values are indexed starting at 0, so results in permutation in S_candidate_length);
	  (*avoidervector[candidate_length + 1]).push_back(extendedperm);
	}
      }
    }

    // Need to add newly found avoiders to next_candidates 
    if (candidate_length + 1 < maxsize) {
      uint64_t tempbits = avoidbits;
      while (tempbits != 0) {
	// i traverses the positions of avoidbits with on-bits from greatest to least
	int i = lastonepos(tempbits); // next position which avoidbits has a one in
	// doing lastonepos first is important for next_candidates ordering
	// so that it's easy to partition next_candidates for recursion
	tempbits = tempbits - (1 << i);

	perm_t extendedinverse =  extendpartialinverse(candidate_inverse, i, candidate_length, maxavoidsize + 1);
	perm_t extendedperm = setdigit(addpos(candidate, i), i, candidate_length); // insert candidate_length in i-th position (remember, values are indexed starting at 0, so results in permutation in S_candidate_length);
	next_candidates.push_back(extendedperm);
	next_candidate_inverses.push_back(extendedinverse);
	if (!justcount) (*avoidervector[candidate_length + 1]).push_back(extendedperm);
	else   *(numavoiders[candidate_length + 1]) += 1;
	bitmaps.push_back(insertbit(avoidbits, i + 1, 1)); // using which insertion positions resulted in an avoider, build avoidbits for each new avoider
      }
    }
  }

  // In the next level of the recursion, we will want to be able to
    // delete any of the largest k letters of the new candidates, and
    // look the result up in the avoidmap. Thus it is okay for us to
    // have fixed the relative order of all but the largest (k-1)
    // letters of the candidates at this level of the
    // recursion. Similarly, we can fix all but the largest (k-1)
    // letters in the next level of the recursion. That means we want
    // to pivot on the position (relative to smaller letters) of the
    // k-th to largest letter of each next_candidate, which is
    // zero-indexed position (candidate_length + 1 - k).  At each
    // level of the recursion, we use n^{k-1} memory There are n
    // levels of the recursion, resulting in n^k total memory
    // usage. In particular, this is a factor of n better than the
    // original memory hack got because we are compressing up to n
    // objects into single bitmaps now.  The recursion partitions the
    // next_candidates based on the pivot.  Note that because of the
    // order we introduce new candidates, currently, all the elements
    // in each part of the partition are already in a row in
    // next_candidates.


  
  if (candidate_length + 1 < maxsize) {
    // If candidate_length + 1 - maxavoidsize \ge 0, we want to do the real
    // recursive step

    if (candidate_length + 1 < maxavoidsize) { 
      buildavoiders_dynamic_helper(candidate_length + 1, patternset, shiftedprefixmap, maxavoidsize, maxsize,  avoidervector, numavoiders, justcount, next_candidates, next_candidate_inverses, 0, next_candidates.size() - 1, avoidermap_write, bitmaps);
    } else {
      // for convenience, we set a variable length = length of common
      // ancestor in tree of candidates of current candidates
      int length = (int)candidate_length - maxavoidsize + 1; 
      // Could be zero if this is the first time we are doing the
      // actual recursion.

      int64_t candidates_prev_end_pos = -1; // index in next_candidates of final candidate in previous part of partition
      int64_t candidates_end_pos = -1; // index of final candidate in current part of partition

      // We are pivoting based on the number of letters to the left of
      // the (length + 1)-th smallest letter that are less than it.  i
      // iterates through the possible values for this pivot. By
      // iterating i backwards, next_candidates is already ordered
      // with based on the pivot value.
      for (int i = length; i >= 0; i--) {
	// cand_pos goes through the indices of the next_candidates with this pivot value
	for (uint64_t cand_pos = candidates_prev_end_pos + 1; cand_pos < next_candidates.size(); cand_pos++) {
	  perm_t next_candidate_inverse = next_candidate_inverses[cand_pos];
	  int lenpos = getdigit(next_candidate_inverse, length);
	  int skip_counter = 0; // counts number of letters larger than length to length's left
	  for (int t = length + 1; t <= candidate_length; t++) {
	    // Note: This factor of k would be a real pain to get rid
	    // of, as far as I can tell (Also there is a factor of k
	    // in the extendedinverse computation).  Consequently,
	    // while the conjecture in the paper would give O(|S_{\le
	    // n - 1}(\Pi)|n) time algorithm (because of the prefix
	    // hack), killing a factor of k from the original
	    // algorithm, the same conjecture seems not to be able to
	    // cut the factor of k off the newer algorithm in which a
	    // factor of n has already been chopped off. Of course,
	    // from a practical standpoint, the prefix hack is still
	    // useful since this for loop has almost nothing inside
	    // it.
	    if (getdigit(next_candidate_inverse, t) < lenpos) skip_counter++;
	  }
	  int norm_counter = lenpos - skip_counter; // position of (length + 1)-th smallest letter relative to smallest length letters of candidate
	  
	  if (norm_counter != i) { // at this point, we've hit the end of the run
	    //displayperm(next_candidate);
	    break; 
	  }
	  candidates_end_pos = cand_pos;
	}
	if (candidates_prev_end_pos + 1 <= candidates_end_pos) {
	  cilk_spawn buildavoiders_dynamic_helper(candidate_length + 1, patternset, shiftedprefixmap,maxavoidsize, maxsize,  avoidervector, numavoiders, justcount, next_candidates, next_candidate_inverses, candidates_prev_end_pos + 1, candidates_end_pos, avoidermap_write, bitmaps);
	  candidates_prev_end_pos = candidates_end_pos;
	}
      }
    }
  }
  cilk_sync;
}


// Dynamic buildavoiders algorithm, implemented using asymptotically memory efficient algorithm, and with Cilk for parallelism
void buildavoiders_dynamic(const hashdb &patternset, int maxavoidsize, int maxsize,  vector < list < perm_t > > &avoidervector, vector < uint64_t > &numavoiders, bool justcount) {
  if (!justcount) avoidervector.resize(maxsize + 1);
  else numavoiders.resize(maxsize + 1);

  hashdb prefixmap(1<<3);
  addprefixes(patternset, prefixmap);

  hashdb cutpatternset(1<< 3); // contains patterns with their largest letter cut out
  hashdb shiftedprefixmap(1<<3); // prefixes for cutpatternset
  vector <perm_t> patterns;
  patternset.getvals(patterns);
  for (int i = 0; i < patterns.size(); i++) {
    perm_t perm = patterns[i];
    int patternlength = getmaxdigit(perm) + 1;
    perm_t inverse = getinverse(perm, patternlength);
    cutpatternset.add(killpos(perm, getdigit(inverse, patternlength - 1)));
  }
  addprefixes(cutpatternset, shiftedprefixmap);


  vector < perm_t > sizetwo;
  vector < perm_t >  sizetwoinverses;
  vector < uint64_t > bitmaps;
  hashmap currentavoiders((1L << 10), 8); // might be worth playing with size of first argument
  
  perm_t startperm = 0;
  uint64_t startpatternmap = 3L;
  currentavoiders.add(startperm, &startpatternmap); // permutations in S_2
  if (!justcount) {
    avoidervector[1].push_back(startperm);
  } else {
    numavoiders[1] = 1;
  }
  perm_t startperm1 = stringtoperm("12");
  perm_t startperm2 = stringtoperm("21");
  if (!justcount) {
    avoidervector[2].push_back(startperm1);
    avoidervector[2].push_back(startperm2); 
  } else {
    numavoiders[2] = 2;
  }
  sizetwo.push_back(startperm1);
  sizetwoinverses.push_back(getinverse(startperm1, 2));
  sizetwo.push_back(startperm2);
  sizetwoinverses.push_back(getinverse(startperm2, 2));
  
  bitmaps.push_back(7L);
  bitmaps.push_back(7L);
  
  vector < cilk::reducer< cilk::op_list_append<perm_t> > > avoidervectortemp(maxsize + 1);
  cilk::reducer< cilk::op_add<uint64_t> > *numavoiderstemp = new cilk::reducer< cilk::op_add<uint64_t> >[maxsize + 1];
  buildavoiders_dynamic_helper(2, patternset, shiftedprefixmap, maxavoidsize, maxsize, avoidervectortemp, numavoiderstemp, justcount, sizetwo, sizetwoinverses, 0, 1, currentavoiders, bitmaps);
  
  for (int i = 3; i < maxsize + 1; i++) {
    if (!justcount) avoidervector[i] = avoidervectortemp[i].get_value();
    // note: we copy the entire avoidervectortemp[i] into avoidervector[i] because
    // it is convenient for the user to not have a data structure which
    // contains a bunch of cilk reducers. However, this is a bit of a waste of work.
    // It's not asymptotically bad though.
    else numavoiders[i] = (uint64_t) numavoiderstemp[i].get_value();
  }
  delete[] numavoiderstemp;
}


// Example:
// string patternlist = "1234 3214"; // space separated list of patterns; need not be same sizes; must be in S_{<10}
// vector < list < perm_t > > avoidervector;
// buildavoidersfrompatternlist(patternlist, 10, avoidervector); // now avoidervector contains S_n(patternlist) stored in avoidervector[n] for 0 < n < 11
void buildavoidersfrompatternlist(string patternlist, int maxpermsize, vector < list < perm_t > > &avoidervector) {
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);
  vector < uint64_t > numavoiders;
  if (USEBRUTE) buildavoiders_brute(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, false, (1L << 10));
  else buildavoiders_dynamic(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, false);
}

// just for use when GETSTAT IS TRUE
uint64_t getstat1() {
  if (!GETSTAT)  return 1;
  return stat1;
}
uint64_t getstat2() {
  if (!GETSTAT)  return 1;
  return stat2;
}
uint64_t getstat3() {
  if (!GETSTAT)  return 1;
  return stat3;
}
uint64_t getstat4() {
  if (!GETSTAT)  return 1;
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
  vector < list < perm_t > > avoidervector;
  if (VERBOSE) cout<<"Effective pattern size "<<maxpatternsize<<endl;
  if (USEBRUTE) buildavoiders_brute(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true, (1L << 10));
  else  buildavoiders_dynamic(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true);
  stat3 = numavoiders[maxpermsize];
  timestamp_t end_time = get_timestamp();
  return (end_time - start_time)/1000000.0;
}

// Example:
// string patternlist = "1234 3214"; // space separated list of patterns; need not be same sizes; must be in S_{<10}
// vector < uint64_t > numavoiders;
// buildavoidersfrompatternlist(patternlist, 10, numavoiders); // now avoidervector contains |S_n(patternlist)| stored in numavoiders[n] for 0 < n < 11.
void countavoidersfrompatternlist(string patternlist, int maxpermsize, vector < uint64_t > &numavoiders) {
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);
  vector < list < perm_t > > avoidervector;
  if (VERBOSE) cout<<"Effective pattern size "<<maxpatternsize<<endl;
  if (USEBRUTE) buildavoiders_brute(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true, (1L << 10));

  else {
    buildavoiders_dynamic(patternset, maxpatternsize, maxpermsize, avoidervector, numavoiders, true);
  }
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
    if (verbose) {
      timestamp_t diff = end_time - start_time;
      double      diff_d = diff* 1e-6;
      cout<< "Time elapsed (s): "<<diff_d<<endl;
    }
    stat3 = numavoiders[maxpermsize];
    if (verbose && GETSTAT && !USEBRUTE) cout<<(double)stat2/(double)stat1<<" prefixes looked at on average per call to isavoider"<<endl;
    if (verbose && GETSTAT && !USEBRUTE) cout<<(double)stat2/(double)stat3<<" prefixes looked at on average per actual avoider"<<endl;
  }
  return;
}



// Parallel version of countavoidersfromfile. (Both use parallel dynamic alg. This one also parallelizes work for different sets of patterns)
void countavoidersfromfile_parallel(ifstream &infile, ofstream &outfile, int maxpermsize, bool verbose) {
  // Note that CLOCK_MONOTONIC isn't defined on Mac OS
  // In other places in the code, we use a different, portable
  // timing mechanism. This is just a quick fix for this part.
  #ifdef CLOCK_MONOTONIC
  struct timespec start_p,end_p;
  clock_gettime(CLOCK_MONOTONIC, &start_p);
  #endif

  string line;
  vector <string> input;
  while (getline(infile, line)) {
    input.push_back(line);
  }
  uint64_t vecsize = input.size();
  vector < vector <uint64_t> > output(vecsize, vector <uint64_t> (0));
  vector < long long > tdiffs(vecsize);
  #ifndef _NO_CILK
  #pragma cilk grainsize = 1
  #endif
  cilk_for (uint64_t ii = vecsize; ii > 0; ii--) {
    uint64_t i = ii-1;
    //cilk_for (uint64_t i = 0; i < vecsize; i++) {
      #ifdef CLOCK_MONOTONIC
    struct timespec start,end;
    vector < uint64_t > numavoiders(maxpermsize + 1);
    countavoidersfrompatternlist(input[i], maxpermsize, numavoiders);
    std::swap(output[i], numavoiders);
    clock_gettime(CLOCK_MONOTONIC, &start);
    clock_gettime(CLOCK_MONOTONIC, &end);
    tdiffs[i] = (end.tv_sec - start.tv_sec)*1000000000 + (end.tv_nsec - start.tv_nsec);
    #endif
    //printf("Another finished \n");
  }
  #ifdef CLOCK_MONOTONIC
  long long maxv = 0, sumv = 0;
  for (auto v : tdiffs) {
    maxv = max(maxv, v);
    sumv += v;
  }
  #endif
  for (uint64_t i = 0; i < vecsize; i++) {
    outfile<<"#"<<input[i]<<endl;
    for (int j = 1; j <= maxpermsize; j++) outfile<<output[i][j]<<" ";
    outfile<<endl;
  }
  #ifdef CLOCK_MONOTONIC
  clock_gettime(CLOCK_MONOTONIC, &end_p);
  long long total_t = (end_p.tv_sec - start_p.tv_sec)*1000000000 + (end_p.tv_nsec - start_p.tv_nsec);
  if (verbose) printf("For total_t=%.6fs maxv=%.6fs sumv=%.6fs parallelism=%f\n", total_t*1e-9, maxv*1e-9, sumv*1e-9, (double)sumv/(double)maxv);
  #endif
  return;
}
