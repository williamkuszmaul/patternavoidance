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
#include "perm.h"
using namespace std;

// for patterns in S_{<10} can use like this:
//  string permlist = "3124 4123 3142 4132";
//  makepatterns(permlist, patternset, maxpatternsize);
void makepatterns(string permlist, hashdb &patternset, int &maxpatternsize) {
  uint64_t lettersused = 0; // Used to detect malformed patterns
  maxpatternsize = 0; 
  int pos = 0;
  perm_t perm = 0;
  //cout<<"Pattern set: "<<endl;
  for (int i = 0; i < permlist.size(); i++) {
    if (permlist[i] == ' ') {
      patternset.add(perm);
      maxpatternsize = max(pos, maxpatternsize);
      if (lettersused + 1 != ((uint64_t)1 << pos)) {
        cout<<"Malformed pattern detected in set "<<permlist<<endl;
        assert (false);
      }
      pos = 0;
      perm = 0;
      lettersused = 0;
    } else {
      lettersused += (1 << (uint64_t)(permlist[i] - '0' - 1));
      perm = setdigit(perm, pos, (int)(permlist[i] - '0' - 1));
      pos++;
    }
  }

  if (permlist[permlist.size() - 1] != ' ') { // if no space at end
    patternset.add(perm);
    if (lettersused + 1 != ((uint64_t)1 << pos)) {
      cout<<"Malformed pattern detected in set "<<permlist<<endl;
      assert (false);
    }
  }
  maxpatternsize = max(pos, maxpatternsize);
}

// adds all the complements of normalizations of prefixes of perm to table
void addprefixeshelper(perm_t perm, int length, hashdb &table) {
  perm_t entry = 0;
  perm_t inverse = getinverse(perm, length);
  uint32_t seenpos = 0; // bit map of which letters we've seen so far
  for (int i = 0; i < length; i++) {
    extendnormalizetop(perm, inverse, length, i, entry, seenpos);
    if (!table.contains(entry)) {
      table.add(entry);
    }
  }
}

// build prefix table containing all complements of normalizations of prefixes of every perm in permset
void addprefixes(const hashdb &permset, hashdb &table) {
  vector <perm_t> patterns;
  permset.getvals(patterns);
  for (int i = 0; i < patterns.size(); i++) {
    perm_t perm = patterns[i];
    int length = getmaxdigit(perm) + 1;
    addprefixeshelper(perm, length, table);
  }
}
