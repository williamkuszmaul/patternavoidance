// THIS FILE IS IN DEVELOPMENT AND IS NOT GUARANTEED TO DO ANYTHING!

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
//#include "wordavoidance.h"
#include "utilities.h"
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

void wordavoiders_raw_dynamic(const hashdb &patternprintset, int maxavoidsize, int maxsize, vector < uint64_t > &numavoiders, int base) {
  hashdb avoidset = hashdb(1000);
  std::queue<perm_t> avoiderstoextend; // queue of avoiders built so far.
  for (int i = 0; i < base; i++)  avoiderstoextend.push(setdigit(0, 0, i));
  for (int i = 0; i < base; i++)  avoidset.add(prepareforavoidset(setdigit(0, 0, i), 1));
  numavoiders[1] = base;

  for (int size = 1; size < maxsize; size++) { // build avoiders of length size + 1
    //cout<<size<<" ----------------"<<endl;
    for (int index = 0; index < numavoiders[size]; index++) {
      perm_t word = avoiderstoextend.front();
      avoiderstoextend.pop();
      for (int letter = 0; letter < base; letter++) {
	perm_t newword = setdigit(word, size, letter);
	bool isavoider = true;

	if (size + 1 <= maxavoidsize) {
	  perm_t wordaspattern = 0;
	  wordaspattern = wordfingerprint(newword, size + 1);
	  if (patternprintset.contains(wordaspattern)) {
	    isavoider = false;
	  }
	}
	if (isavoider) {
	  for (int killedpos = 0; killedpos < size + 1 && killedpos <= maxavoidsize; killedpos++) {
	    perm_t shrunkword = killpos(newword, killedpos);
	    if (!avoidset.contains(prepareforavoidset(shrunkword, size))) {
	      isavoider = false;
	      break;
	    }
	  }
	}
	if (isavoider) {
	  avoiderstoextend.push(newword);
	  avoidset.add(prepareforavoidset(newword, size + 1));
	  // cout<<size + 1<<endl;
	  //displayperm(newword);
	  numavoiders[size + 1]++;
	}
      }
    }
  }
}

// REQUIREMENTS:
// patternlist contains patterns of size 8 or less (b/c of fingerprint function)
// base <= 31 and base <= max digit size in word (first requirement because of fingerprint function)
// Cannot have word comprised of 16 16s if using 64 bit permutation representation (hash table -1 error)
// maxwordsize <= number of digits we can store in word
void countavoidersfrompatternlist(string patternlist, int base, int maxwordsize, vector < uint64_t > &numavoiders) {
  hashdb patternprintset(1);
  numavoiders.resize(maxwordsize + 1);
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

  wordavoiders_raw_dynamic(patternprintset, maxpatternsize, maxwordsize, numavoiders, base);
}

int main() {
  int maxwordsize = 10;
  int base = 7;
  string patternlist = "1234";
  vector < uint64_t > numavoiders; //(11);
  countavoidersfrompatternlist(patternlist, base, maxwordsize, numavoiders);
  //wordavoiders_raw_dynamic(patternprintset, 2, 10, numavoiders, 2);
  for (int i = 1; i <= maxwordsize; i++) {
    cout<<numavoiders[i]<<" ";
  }
  cout<<endl;
  return 0;
}
