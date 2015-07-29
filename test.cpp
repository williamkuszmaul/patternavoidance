#include <assert.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <time.h> 
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <stdint.h>
// #include <unordered_set>
#include <queue>
#include "hashdb.h"
#include <sys/time.h>
#include "fastavoidance.h"
using namespace std;


// for patterns in S_{<10} can use like this:
//  string permlist = "3124 4123 3142 4132";
//  makepatterns(permlist, patternset);
uint64_t makepatterns(string permlist, hashdb &patternset) {
  int pos = 0;
  uint64_t perm = 0;
  cout<<"Pattern set: "<<endl;
  for (int i = 0; i < permlist.size(); i++) {
    if (permlist[i] == ' ') {
      patternset.add(perm);
      displayperm(perm);
      pos = 0;
      perm = 0;
    } else {
      perm = setdigit(perm, pos, (int)(permlist[i] - '0' - 1));
      pos++;
    }
  }
  // No space at end
  patternset.add(perm);
  displayperm(perm);
      
  cout<<endl;
}

int main() {
  int maxpatternsize = 4;
  int permsize = 8;
  assert(permsize <= 16);
  hashdb patternset = hashdb(1<<3);
  cout<<"Avoid set: ";

  string permlist = "3124 4123 3142 4132";
  makepatterns(permlist, patternset);
  // Note: could have boolian for whether to fill avoiders vector. And could have boolian for whether to count avoiders of each size. Or just a separte function to count avoiders of each size.

  //perm = setdigit(perm, 3, 3);
  timestamp_t start_time = get_timestamp();
  vector < int > numavoiders;
  countavoiders(patternset, maxpatternsize, permsize, numavoiders);
  for (int i = 2; i <= permsize; i++) {
    cout<<"Number of avoiders of size "<<i<<" is "<<numavoiders[i]<<endl;
  }
  //  assert(numavoiders == avoidersvector[permsize].size());
  timestamp_t end_time = get_timestamp();
  cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  return 0;
}
