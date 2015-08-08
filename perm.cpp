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
//  makepatterns(permlist, patternset);
uint64_t makepatterns(string permlist, hashdb &patternset, int &maxpatternsize) {
  maxpatternsize = 0; // IN COUNT AVOIDERS AM I CORRECTLY INITIALIZING ARRAY THINGS AT ZERO
  int pos = 0;
  uint64_t perm = 0;
  //cout<<"Pattern set: "<<endl;
  for (int i = 0; i < permlist.size(); i++) {
    if (permlist[i] == ' ') {
      patternset.add(perm);
      //displayperm(perm);
      maxpatternsize = max(pos + 1, maxpatternsize);
      pos = 0;
      perm = 0;
    } else {
      perm = setdigit(perm, pos, (int)(permlist[i] - '0' - 1));
      pos++;
    }
  }

  if (permlist[permlist.size() - 1] != ' ') { // if no space at end
    patternset.add(perm);
  }
  maxpatternsize = max(pos + 1, maxpatternsize);
}



unsigned long long permtonum(uint64_t perm, int length) {
  int answer = 0;
  unsigned long long fac = 1;
  uint64_t seenletters = 0; // bit map of which letters we've seen so far
  for (int i = length - 1; i >= 0; i--) {
    unsigned long long facdig = 0;
    int digit = getdigit(perm, i);
    if (digit != 0) {
      unsigned int temp = seenletters << (32 - digit);
      facdig = __builtin_popcount(temp);
    }
    answer += facdig * fac;
    fac *= (length - i);
    seenletters = seenletters | (1 << digit);
  }
  return answer;
}
