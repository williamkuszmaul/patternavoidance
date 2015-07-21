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


int main() {
    int maxpatternsize = 4;
  int permsize = 12;
  assert(permsize <= 16);
  uint64_t perm = 0;
  perm = setdigit(perm, 0, 0);
  perm = setdigit(perm, 1, 2);
  perm = setdigit(perm, 2, 1);
  perm = setdigit(perm, 3, 3);
  hashdb patternset = hashdb(1<<3);
  patternset.add(perm);
  timestamp_t start_time = get_timestamp();
  cout<<"Avoid set: ";
  displayperm(perm);
  vector < vector < uint64_t > > avoidersvector;
  uint64_t numavoiders = buildavoiders(patternset, maxpatternsize, permsize, avoidersvector);
  cout<<"Number of avoiders of size "<<permsize<<" is "<<numavoiders<<endl;
  assert(numavoiders == avoidersvector[permsize].size());
  timestamp_t end_time = get_timestamp();
  cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  return 0;
}
