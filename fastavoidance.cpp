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
using namespace std;


void displayperm(uint64_t perm) {
  for (int i = 0; i < 16; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

void displayperm(uint64_t perm, int size) {
  for (int i = 0; i < size; i++) cout<<getdigit(perm, i)<<" ";
  cout<<endl;
}

static bool isavoider(uint64_t perm, int maxavoidsize, int length, const hashdb &avoidset, const hashdb &patternset) { 
  uint64_t inverse = getinverse(perm, length);
  if (length <= maxavoidsize && patternset.contains(perm)) { // if is in set of bad patterns
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
      //displayperm(currentperm)
      // could make below req include i != length - 1 && 
      if (!avoidset.contains(currentperm)) {
	return false; // found a subword not avoiding the patterns
      }
    }
  }
  return true;
}


void buildavoiders(const hashdb &patternset, int maxavoidsize, int maxsize,  vector < vector < uint64_t > > &avoidervector, uint64_t plannedavoidsetsize) {
  avoidervector.resize(maxsize + 1);
  
  hashdb avoidset = hashdb(plannedavoidsetsize);
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
    for (int i = 0; i < currentlength + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentlength);
      if (isavoider(extendedperm, maxavoidsize, currentlength + 1, avoidset, patternset)) {
	avoidervector[currentlength + 1].push_back(extendedperm);
      	if (currentlength + 1 < maxsize) {
	  avoiderstoextend.push(extendedperm);
	  avoidset.add(extendedperm);
	  numnextlength++;
	}
      }
    }
  }
}



void countavoiders(const hashdb &patternset, int maxavoidsize, int maxsize, vector < uint64_t > &numavoiders, uint64_t plannedavoidsetsize) {
  numavoiders.resize(maxsize + 1); // counts number of avoiders of size i
  
  hashdb avoidset = hashdb(plannedavoidsetsize);
  uint64_t startperm = 0;
  avoidset.add(startperm); // identity in S_1
  
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
    for (int i = 0; i < currentlength + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentlength);
      if (isavoider(extendedperm, maxavoidsize, currentlength + 1, avoidset, patternset)) {
	numavoiders[currentlength + 1]++;
      	if (currentlength + 1 < maxsize) {
	  avoiderstoextend.push(extendedperm);
	  avoidset.add(extendedperm);
	  numnextlength++;
	}
      }
    }
  }
}


// for patterns in S_{<10} can use like this:
//  string permlist = "3124 4123 3142 4132";
//  makepatterns(permlist, patternset);
uint64_t makepatterns(string permlist, hashdb &patternset, int &maxpatternsize) {
  maxpatternsize = 0;
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


void buildavoidersfrompatternlist(string patternlist, int maxpermsize, vector < vector < uint64_t > > &avoidervector) {
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);
  buildavoiders(patternset, maxpatternsize, maxpermsize, avoidervector, (1L << 10)); // for large cases, make last argument much larger!
}


void countavoidersfrompatternlist(string patternlist, int maxpermsize, vector < uint64_t > &numavoiders) {
  int maxpatternsize;
  hashdb patternset = hashdb(1<<3);
  makepatterns(patternlist, patternset, maxpatternsize);
  countavoiders(patternset, maxpatternsize, maxpermsize, numavoiders, (1L << 10)); // for large cases, make last argument much larger!
}

void countavoidersfromfile(ifstream &infile, ofstream &outfile, int maxpermsize, bool verbose) {
  string line;
  while (getline(infile, line)) {
    outfile<<"#"<<line<<endl;
    vector < uint64_t > numavoiders;
    timestamp_t start_time = get_timestamp();
    countavoidersfrompatternlist(line, maxpermsize, numavoiders);
    timestamp_t end_time = get_timestamp();
    if (verbose) cout<<line<<endl;
    if (verbose) cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
    for (int i = 1; i < numavoiders.size(); i++) {
      outfile<<numavoiders[i]<<" ";
    }
    outfile<<endl;
  }
  return;
}

// ------------------- A different version of countavoiders which uses a bitmap to get some speed-up, based on a hack in permlab

// for special usage in countavoidersv2
static bool isavoiderv2(uint64_t perm, int maxavoidsize, int length, const hashdb &avoidset, const hashdb &patternset) { 
  uint64_t inverse = getinverse(perm, length);
  if (length <= maxavoidsize && patternset.contains(perm)) { // if is in set of bad patterns
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
      //displayperm(currentperm)
      // could make below req include i != length - 1 && 
      if (i < length - 2 && !avoidset.contains(currentperm)) {
	return false; // found a subword not avoiding the patterns
      }
    }
  }
  return true;
}

inline uint64_t getbit(uint64_t word, int pos) {
  return (word >> pos) & 1;
}

inline uint64_t setbit(uint64_t word, uint64_t pos, uint64_t val) { // val is 1 or 0
  return (word & (~(0L) - (1<<pos))) + (1<<pos)*val;
}

// Note that due to bit-shift setup, is not capable to inserting in final position of word
inline uint64_t insertbit(uint64_t word, uint64_t pos, uint64_t val) { // val is 1 or 0
  return (word & ((1 << pos) - 1) ) + ((word >> pos) << (pos + 1)) + val * (1 << pos);
}

// Incorporates the bitmap idea from permlab
// still experimental
void countavoidersv2(const hashdb &patternset, int maxavoidsize, int maxsize, vector < uint64_t > &numavoiders, uint64_t plannedavoidsetsize) {
  numavoiders.resize(maxsize + 1); // counts number of avoiders of size i
  
  hashdb avoidset = hashdb(plannedavoidsetsize);
  uint64_t startperm = 0;
  avoidset.add(startperm); // identity in S_1
  
  std::queue<unsigned long long> avoiderstoextend;
  std::queue<unsigned long long> bitmaps;
  avoiderstoextend.push(startperm);
  bitmaps.push(3L);

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
    bitmaps.pop();
    numleftcurrentlength--;
    for (int i = 0; i < currentlength + 1; i++) {
      if (getbit(bitmap, i) == 1) {
	uint64_t extendedperm = setdigit(addpos(perm, i), i, currentlength);
	if (isavoiderv2(extendedperm, maxavoidsize, currentlength + 1, avoidset, patternset)) {
	  numavoiders[currentlength + 1]++;
	  if (currentlength + 1 < maxsize) {
	    avoiderstoextend.push(extendedperm);
	    avoidset.add(extendedperm);
	    numnextlength++; 
	  }
	} else {
	  bitmap = setbit(bitmap, i, 0);
	}
      }
    }
    if (currentlength + 1 < maxsize) {
      for (int i = 0; i < currentlength + 1; i++) {
	if (getbit(bitmap, i) == 1) {
	  bitmaps.push(insertbit(bitmap, i + 1, 1));
	}
      }
    }
  }
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
					   
