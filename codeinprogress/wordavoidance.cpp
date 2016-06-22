// THIS FILE IS IN DEVELOPMENT AND IS NOT GUARANTEED TO DO ANYTHING!
// I'm using it to run some experiments.

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
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <cilk/reducer_list.h>

#include <queue>
#include <sys/time.h>
#include <vector>
#include "hashdb.h"
#include "hashmap.h"
//#include "wordavoidance.h"
#include "permutilities.h"
#include "perm.h"
using namespace std;


// THINGS TO BE CAREFUL ABOUT:
// Words could now take val -1 if they are large, screwing up hash table. That includes fingerprint (but is not an issue in our current implementation for fingerprint).
// to distinguish words, DO need to know size. Thus avoidset and word fingerprint BOTH include this in their representatinos


// Returns the unique fingerprint of pattern represented by word
// Requires
// letters of value <= 31
// More importantly, never more than 7 letters preceeding another letter and larger than it (assuming uint64_t perm)
// Also size of word  + 1 < max allowed size of word
// SUFFICES to have letters of value <= 31 and word size 8 or less
perm_t wordfingerprint(perm_t word, int size) {
  perm_t print = setdigit(0, size, 1); // so that finger print marks size
  uint32_t seenletters = 0;
  for (int i = 0; i < size; i++) {
    unsigned long long reldig = 0;
    int digit = getdigit(word, i);
    unsigned int temp = seenletters >> (digit);
    reldig = 2*__builtin_popcount(temp);
    if ((seenletters & (1 << digit))  != 0) reldig--;
    // now relcount is count 2 * number larger digits + 1 if current digit also appears
    // have to be careful to do this because can have letters appearing multiple times
    seenletters = seenletters | (1 << digit);
    print = setdigit(print, i, reldig);
  }
  return print;
}

perm_t prepareforavoidset(perm_t word, int length) {
  return setdigit(word, length, 1); // mark length
}

void wordavoiders_raw_dynamic_helper(const hashdb &patternprintset,  int maxavoidsize, int maxsize,  cilk::reducer< cilk::op_add<uint64_t> > *numavoiders, int base, int candidate_length, const vector <perm_t> &candidates, int first_cand_index, int last_cand_index, const hashmap &avoidmap_read) {
  hashmap avoidmap_write(1<<8, 8);
  uint64_t trivialmap = (((uint64_t)1) << base) - 1;
  vector <perm_t> avoiderstoextend;
  for (int candidate_index= first_cand_index; candidate_index <= last_cand_index; candidate_index++) {
    perm_t word = candidates[candidate_index];
    uint64_t avoidbits = trivialmap; // first base bits are ones
    for (int killedpos = candidate_length - 1; killedpos > candidate_length - 1 - maxavoidsize && killedpos >= 0; killedpos--) {
      perm_t shrunkword = killpos(word, killedpos);
      uint64_t bitmap_temp = *((uint64_t *)(avoidmap_read.getpayload(prepareforavoidset(shrunkword, candidate_length - 1))));
      avoidbits = avoidbits & bitmap_temp;
    }
    // still need to check if each of detected avoider might actually be a pattern:
    if (candidate_length + 1 <= maxavoidsize) {
      // Not worth messing with bithacks to make this loop more efficient
      // since we only enter this if statement in small cases.
      for (int addval = 0; addval < base; addval++) {
	if ((avoidbits & (1 << addval)) != 0) {
	  perm_t newword = setdigit(word, candidate_length, addval);
	  perm_t wordaspattern = 0;
	  wordaspattern = wordfingerprint(newword, candidate_length + 1);
	  if (patternprintset.contains(wordaspattern)) {
	    avoidbits = avoidbits - (1 << addval); // was not an avoider after all
	  } 
	}
      }
    }
    if (candidate_length + 1 < maxsize) avoidmap_write.add(prepareforavoidset(word, candidate_length), &avoidbits);
    *(numavoiders[candidate_length + 1]) += __builtin_popcount((uint32_t)avoidbits);
    if (candidate_length + 1 < maxsize) {
      for (int addval = 0; addval < base; addval++) {
	// Probably not worth messing around with bit hacks to improve asymptoics on this for loop.
	// In particular, if number of avoiders grows even sort of fast, then for the
	// avoiders we push onto avoiderstoextend we're going to be doing enough work
	// to amortize the wasted cost in this for loop.
	if ((avoidbits & (1 << addval)) != 0) {
	  perm_t newword = setdigit(word, candidate_length, addval);
	  avoiderstoextend.push_back(newword);
	  //avoidset.add(prepareforavoidset(newword, candidate_length + 1));
	  //cout<<candidate_length + 1<<endl;
	  //displayperm(newword);
	}
      }
    }
  }

  if (candidate_length + 1 < maxsize) {
    // In the next level of the recursion, we will want to be able to
    // delete any of the final k letters of the new candidates, and
    // look the result up in the avoidmap. Thus it is okay for us to
    // have fixed all but the final (k-1) letters of the candidates at
    // this level of the recursion. Similarly, we can fix all but the
    // final (k-1) letters in the next level of the recursion. That
    // means we want to pivot on the k-th to final letter of each
    // next-candidate, which is zero-indexed position
    // (candidate_length + 1 - k).
    // At each level of the recursion, we use n^{k-1} memory
    // There are n levels of the recursion, resulting in n^k
    // total memory usage. In particular, this is a factor of n
    // better than the original memory hack got because we are
    // compressing up to n objects into single bitmaps now.

    // If candidate_length + 1 - maxavoidsize \ge 0, we want to do the real
    // recursive step
    if (candidate_length < maxavoidsize - 1) {
      wordavoiders_raw_dynamic_helper(patternprintset, maxavoidsize, maxsize, numavoiders, base, candidate_length + 1, avoiderstoextend, 0, avoiderstoextend.size() - 1, avoidmap_write);
    } else {
      int next_candidate_index = 0;
      for (int splitval = 0; splitval < base; splitval++) {
	int start_index = next_candidate_index;
	int end_index = next_candidate_index - 1;

	// We use the letter in zero-indexed position candidate_length + 1 - k as our pivot 
	while (next_candidate_index < avoiderstoextend.size()
	       && getdigit(avoiderstoextend[next_candidate_index], candidate_length + 1 - maxavoidsize) == splitval) {
	  end_index++;
	  next_candidate_index++;
	}
	cilk_spawn wordavoiders_raw_dynamic_helper(patternprintset, maxavoidsize, maxsize, numavoiders, base, candidate_length + 1, avoiderstoextend, start_index, end_index, avoidmap_write);
      }
      cilk_sync;
    }
  }
}


// PATTERNS HAVE TO BE AT LEAST SIZE THREE DUE TO HOW AVOIDMAP IS SET UP
void wordavoiders_raw_dynamic(const hashdb &patternprintset, int maxavoidsize, int maxsize,  cilk::reducer< cilk::op_add<uint64_t> > *numavoiders, int base) {
  hashmap avoidmap(1<<8, 8);
  std::vector<perm_t> avoiderstoextend; // queue of avoiders built so far.
  //for (int i = 0; i < base; i++)  avoidset.add(prepareforavoidset(setdigit(0, 0, i), 1));
  uint64_t trivialmap = (((uint64_t)1) << base) - 1;
  for (int i = 0; i < base; i++)  avoidmap.add(prepareforavoidset(setdigit(0, 0, i), 1), &trivialmap);
  *(numavoiders[1]) += base;
  for (int i = 0; i < base; i++) {
    for (int j = 0; j < base; j++) {
      perm_t temp_perm = setdigit(setdigit(0, 0, i), 1, j);
      avoiderstoextend.push_back(prepareforavoidset(temp_perm, 2));
      //avoidset.add(prepareforavoidset(temp_perm, 2));
      avoidmap.add(prepareforavoidset(temp_perm, 2), &trivialmap);
    }
  }
  *(numavoiders[2]) += base * base;

  wordavoiders_raw_dynamic_helper(patternprintset, maxavoidsize, maxsize, numavoiders, base,  2, avoiderstoextend, 0, avoiderstoextend.size() - 1, avoidmap);
}

// REQUIREMENTS:
// patternlist contains patterns of size 8 or less (b/c of fingerprint function)
// base <= 31 and base <= max digit size in word (first requirement because of fingerprint function)
// Cannot have word comprised of 16 16s if using 64 bit permutation representation (hash table -1 error)
// maxwordsize <= number of digits we can store in word
void countavoidersfrompatternlist(string patternlist, int base, int maxwordsize, vector <uint64_t> &numavoiders) {
  hashdb patternprintset(1);
  numavoiders.resize(maxwordsize + 1);
  cilk::reducer< cilk::op_add<uint64_t> > numavoiderstemp[maxwordsize + 1];
  int maxpatternsize = 0;
  int pos = 0;
  perm_t word = 0;
  //cout<<"Pattern set: "<<endl;
  for (int i = 0; i < patternlist.size(); i++) {
    if (patternlist[i] == ' ') {
      patternprintset.add(wordfingerprint(word, pos));
      maxpatternsize = max(pos, maxpatternsize);
      pos = 0;
      word = 0;
    } else {
      word = setdigit(word, pos, (int)(patternlist[i] - '0' - 1));
      pos++;
    }
  }
  if (patternlist[patternlist.size() - 1] != ' ') { // if no space at end
    patternprintset.add(wordfingerprint(word, pos));
  }
  maxpatternsize = max(pos, maxpatternsize);

  wordavoiders_raw_dynamic(patternprintset, maxpatternsize, maxwordsize, numavoiderstemp, base);
  for (int i = 0; i < maxwordsize + 1; i++) {
    numavoiders[i] = numavoiderstemp[i].get_value();
  }
}

class Bitmap {
public:
  vector <bool> map;
  int size;
  Bitmap(int requestedsize) {
    size = requestedsize;
    map.resize(size);
  }
  void setpos(int pos) {
    map[pos] = true;
  }
  void clearpos(int pos) {
    map[pos] = false;
  }
  bool getpos(int pos) {
    return map[pos];
  }
  void display() {
    for (int i = 0; i < size; i++) cout<<map[i]<<" ";
    cout<<endl;
  }
  int numonbits(){
    int answer = 0;
    for (int i = 0; i < size; i++) answer += (int)map[i];
    return answer;
  }
  void clear() {
    for (int i = 0; i < size; i++) map[i] = false;
  }
};


// replaces map with next bitmap lexicographically to have same number of on bits as this one
// returns false if map is the largest such map
static bool getnextmap(Bitmap &map) {
  int i = map.size - 1;
  int numleadingones = 0;
  while (i >= 0 && map.getpos(i) == 1) {
    numleadingones++;
    map.clearpos(i);
    i--;
  }
  while (i >= 0 && map.getpos(i) == 0) {
    i--;
  }
  if (i < 0) return false; 
  map.clearpos(i);
  i++;
  for (int j = 0; j < numleadingones + 1; j++) {
    map.setpos(i+j);
  }
  return true;
}

void allsubsetsofsize(string patterns[], int totalsize,  int setsize, int maxwordsize, int base, string outfile) {
  ofstream outstream;
  outstream.open(outfile, std::ofstream::trunc);
  vector < string > patternsets;
  Bitmap map(totalsize);
  for (int j = 0; j < setsize; j++) map.setpos(j);

  // First build sets of patterns
  while (1)  {
    string patternset = "";
    for (int i = 0; i < totalsize; i++) {
      if (map.getpos(i)) patternset = patternset + " " + patterns[i];
    }
    patternsets.push_back(patternset);
    if (!getnextmap(map)) break;
  }

  // Next compute number of avoiders for each set of patterns
  vector < vector <uint64_t> > sequences(patternsets.size(), vector <uint64_t> (maxwordsize + 1));
  //cout<<patternsets.size()<<endl;
#pragma cilk grainsize = 1
  cilk_for (int i = 0; i < patternsets.size(); i++) {
    vector <uint64_t> numavoiders(maxwordsize + 1);
    //cout<<pattersets[i]<<endl;
    countavoidersfrompatternlist(patternsets[i], base, maxwordsize, numavoiders);
    for (int x = 1; x <= maxwordsize; x++) {
      sequences[i][x] = numavoiders[x];
    }
  }

  // Finally, report results
  for (int i = 0; i < patternsets.size(); i++) {
    outstream<<patternsets[i]<<endl;
    for (int x = 1; x <= maxwordsize; x++) {
      outstream<<sequences[i][x]<<" ";
    }
    outstream<<endl;
  }
  outstream.close();
}


int main() {
  string array[] = {"1111","1112","1121","1122","1123","1211","1212","1213","1221","1222","1223","1231","1232","1233","1243","1312","1321","1322","1323","1324","1332","1342","1423","1432","2111","2112","2113","2121","2122","2123","2131","2132","2133","2134","2143","2211","2212","2213","2221","2231","2311","2312","2313","2314","2321","2331","2341","2413","2431","3112","3121","3122","3123","3124","3132","3142","3211","3212","3213","3214","3221","3231","3241","3312","3321","3412","3421","4123","4132","4213","4231","4312","4321"};
  int base = 7;
  int maxsize = 8;
  int minsetsize = 4;
  int maxsetsize = 4;
  for (int i = minsetsize; i <= maxsetsize; i++) {
    cout<<"Starting sets of size "<<i<<endl;
    allsubsetsofsize(array, sizeof(array) / sizeof(char**), i, maxsize, base, "worddata" + to_string(i));
  }
  cout<<"Done"<<endl;
  return 0;
}
