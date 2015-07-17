// simplest make is g++ -O2 -std=c++11  fastavoidance.cpp  -g -o fastavoidance
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

#define ULL (uint64_t)
#define one ((uint64_t) 1)
#define zero ((uint64_t) 0)
#define bits4 ((uint64_t) 15)

void displayperm(uint64_t perm) {
  for (int i = 0; i < 16; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

void displayperm(uint64_t perm, int size) {
  for (int i = 0; i < size; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

bool isavoider(uint64_t perm, int maxavoidsize, int length, hashdb &avoidset, hashdb &patternset) { 
  uint64_t inverse = getinverse(perm, length);
  if (patternset.contains(perm)) { // if is in set of bad patterns
      return false;
  }
  uint64_t currentperm = perm;
  if (length > 1) { // don't deal with permutations of size zero
    for (int i = length - 1; i >= 0 && i >= length - maxavoidsize - 1; i--) {
      if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
        currentperm = addpos(currentperm, getdigit(inverse, i + 1));
	currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
      }
      currentperm = killpos(currentperm, getdigit(inverse, i));
      //displayperm(currentperm);
      if (!avoidset.contains(currentperm)) {
	return false; // found a subword not avoiding the patterns
      }
    }
  }
  return true;
}

long long buildavoiders(hashdb &patternset, int maxavoidsize, int maxsize, vector < vector < uint64_t > > &avoidervector) {
  int numavoiders = 0;
  avoidervector.resize(maxsize + 1); // avoidervector[i] will contain the avoiders of size i. [avoidervector[0] will be empty]
  
  hashdb avoidset = hashdb(1<<20);
  uint64_t startperm = 0;
  avoidset.add(startperm); // identity in S_1
  avoidervector[1].push_back(startperm);
  
  std::queue<long long> avoiderstoextend;
  avoiderstoextend.push(startperm);

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
    avoiderstoextend.pop();
    numleftcurrentlength--;
    if (currentlength >= maxsize) {
      break;
    }    
    for (int i = 0; i < currentlength + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentlength);
      if (isavoider(extendedperm, maxavoidsize, currentlength + 1, avoidset, patternset)) {
      	if (currentlength + 1 == maxsize) {
	  numavoiders++;
	}
	avoidervector[currentlength + 1].push_back(extendedperm);
      	avoiderstoextend.push(extendedperm);
	numnextlength++;
	avoidset.add(extendedperm);
      }
    }
  }
  return numavoiders;
}




// int main() {  
//   int maxpatternsize = 3;
//   int permsize = 14;
//   assert(permsize <= 16);
//   uint64_t perm = 0;
//   perm = setdigit(perm, 0, 0);
//   perm = setdigit(perm, 1, 2);
//   perm = setdigit(perm, 2, 1);
//   //perm = setdigit(perm, 3, 3);
//   hashdb patternset = hashdb(1<<3);
//   patternset.add(perm);
//   timestamp_t start_time = get_timestamp();
//   cout<<"Avoid set: ";
//   displayperm(perm);
//   vector < vector < uint64_t > > avoidersvector;
//   uint64_t numavoiders = buildavoiders(patternset, maxpatternsize, permsize, avoidersvector);
//   cout<<"Number of avoiders of size "<<permsize<<" is "<<numavoiders<<endl;
//   assert(numavoiders == avoidersvector[permsize].size());
//   timestamp_t end_time = get_timestamp();
//   cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
//   return 0;
// }
					   
