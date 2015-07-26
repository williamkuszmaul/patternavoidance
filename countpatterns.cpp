// g++ -O0 -std=c++11 -g countpatterns.cpp hashdb.h hashdb.cpp fastavoidance.h fastavoidance.cpp -o countpatterns
// Note: Valgrind falsely reports 72,704 bytes in use at exit, even when I return at beginning of main. This seems to be a documented valgrind bug.

// scp countpatterns ec2:./
// ssh ec2
// /usr/bin/time ./countpatterns

// top gets things like gb usage. run in different shell
// can get big speedup by planning hash table size appropriately

//ssh ec2 to get into ec2
//scp countpatterns ec2:./
#include <assert.h>
#include <string.h>
//#include <iostream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <time.h> 
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <stdint.h>
#include <unordered_map>
#include <queue>
#include "hashdb.h"
#include "hashmap.h"
#include <sys/time.h>
#include "fastavoidance.h"
using namespace std;

const int permsize = 12;
const int maxpatternsize = 5;

vector < vector<int> > finaldata;

void initfinaldata() {
  finaldata.resize(permsize + 1);
  int fac = 1;
  for (int i = 1; i <= permsize; i++) {
    fac *= i;
    finaldata[i].resize(fac);
  }
  return;
}

template <unsigned int MAXPATTERNSIZE> 
class myarray {
  public:
    unsigned short dataset[MAXPATTERNSIZE];
    myarray() {}
    myarray(unsigned short *array) {
    for (int i = 0 ; i < MAXPATTERNSIZE; i++) {
      dataset[i] = array[i];
    }
  }
};


// This is a non-portable way to use the built-in popcnt. But using -march=native in makefile and builtin g++ command seems to work just as well and be more portable
// inline int popcount (unsigned int input) {
//   int answer;
//   asm("popcnt %1, %0" : "=r" (answer) : "r" (input));
//   return answer;
// }

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



inline void setPvals(unsigned long long perm, unsigned short* payload, hashmap &Phashmap) {
  Phashmap.add(perm, payload);
}

inline unsigned short getPval(unsigned long long perm, int i, hashmap &Phashmap) {
  return ((unsigned short*)Phashmap.getpayload(perm))[i];
}

static void Pcount(uint64_t perm, int length, const hashdb &patternset, hashmap &Phashmap) { 
  uint64_t inverse = getinverse(perm, length);
  unsigned short Pvals[maxpatternsize + 2]; // vals range from [0...maxpatternsize+1]
  int Pvalpos = maxpatternsize + 1;
  while (Pvalpos > maxpatternsize || Pvalpos > length) {
    Pvals[Pvalpos] = 0;
    Pvalpos--;
  }
  if (Pvalpos == length) {
    if (patternset.contains(perm)) Pvals[Pvalpos] = 1;
    else Pvals[Pvalpos] = 0;
    Pvalpos--;
  }

  // Now Pvalpos < length and <= maxpatternsize. So the recurrence starts to apply

  unsigned short oldPvals[Pvalpos + 1];
  int oldPvalpos = 0;
  uint64_t currentperm = perm;
  assert(length > 1); // don't want to deal with length 1 case
  for (int i = length - 1; i >= 0 && i >= length - maxpatternsize - 1; i--) {
    if (i < length - 1) { // add back in digit we deleted a moment ago, but with value one smaller
      currentperm = addpos(currentperm, getdigit(inverse, i + 1));
      currentperm = setdigit(currentperm, getdigit(inverse, i + 1), i);
    }
    currentperm = killpos(currentperm, getdigit(inverse, i));
    oldPvals[oldPvalpos] = getPval(currentperm, length - 1 - i, Phashmap);
    //assert(getPval(currentperm, length - 1 - i, Pmap) ==
    //   getPval2(currentperm, length - 1- i, Phashmap));
    oldPvalpos++;
  }

  while (Pvalpos >= 0) {
    Pvals[Pvalpos] = Pvals[Pvalpos + 1] + oldPvals[Pvalpos];
    Pvalpos--;
  }
    if (length != permsize) {
      //setPvals(perm, Pvals, Pmap);
      setPvals(perm, Pvals, Phashmap);
    }
    int index = permtonum(perm, length);
    finaldata[length][index] = Pvals[0];
  //  displayperm(perm);
  //cout<<"num hits: "<<Pvals[0]<<endl;
}


void buildpermutations(uint64_t perm, int currentsize, int finalsize, hashdb &patternset, hashmap &Phashmap) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize, patternset, Phashmap);
    }
  } else {
    Pcount(perm, currentsize, patternset, Phashmap);
  }
}

// LENGTH IS FUNNY SOME TIMES!
void createPmap(uint64_t finalsize, hashdb &patternset, timestamp_t start_time, hashmap &Phashmap) {
  unsigned short temp[4] = {0, 0, 0, 0};
  //setPvals(0, temp, Pmap);
  setPvals(0, temp, Phashmap);
  uint64_t perm1 = setdigit(0L, 0, 1L);
  uint64_t perm2 = setdigit(0L, 1, 1L);
  for (int i = 2; i <= finalsize; i++) {
    buildpermutations(0L, 1, i, patternset, Phashmap);
    timestamp_t current_time = get_timestamp();
    cout<< "Time elapsed to build perms of size "<<i<<" in seconds: "<<(current_time - start_time)/1000000.0L<<endl;
  }
  return;
}

unsigned long long factorial(long long num) {
  int answer = 1;
  while (num > 1) {
    answer *= num;
    num--;
  }
  return answer;
}

int main() {
  initfinaldata();
  //  void initpopCountOfByte256();
  assert(permsize <= 16);
  uint64_t perm = 0;
  perm = setdigit(perm, 0, 0);
  perm = setdigit(perm, 1, 1);
  perm = setdigit(perm, 2, 2);
  perm = setdigit(perm, 3, 3);
  perm = setdigit(perm, 4, 4);
  hashdb patternset = hashdb(1<<3);
  hashdb patternset2 = hashdb(1<<3);
  patternset.add(perm);
  patternset2.add(perm);
 
  unsigned long long reservedspace = 0;
  for (int i = 1; i <= permsize - 1; i++) reservedspace += factorial(i);
  hashmap Phashmap(reservedspace * 3, sizeof(myarray<maxpatternsize + 1>));
  timestamp_t start_time = get_timestamp();
  cout<<"Pattern set: ";
  displayperm(perm);
  createPmap(permsize, patternset, start_time, Phashmap);
  timestamp_t end_time = get_timestamp();
  cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  return 0;
}
