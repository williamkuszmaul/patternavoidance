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
#include "perm.h"
using namespace std;




// If the i-th to last letter (starting with i = 1) is greater than j letters to its right, then we add j * (i-1)! to our final answer.
// This yields a unique number because of how base factorial works (https://en.wikipedia.org/wiki/Factorial_number_system)
unsigned long long permtonum(perm_t perm, int length) {
  int answer = 0;
  unsigned long long fac = 1; // will be updated to take values of factorials
  uint64_t seenletters = 0; // bit map of which letters we've seen so far
  for (int i = length - 1; i >= 0; i--) {
    unsigned long long facdig = 0; // will be the number of letters to the right of the i-th letter with value less than the i-th letter
    int digit = getdigit(perm, i);
    if (digit != 0) { 
      unsigned int temp = seenletters << (32 - digit); // Note: shifting by 32 is ill-defined, which is why we explicitly eliminate digit = 0 case.
      facdig = __builtin_popcount(temp);  // By taking the popcount of temp, we get the number of letters to the right of the i-th letter with value less than the i-th letter
      // Note: builtin_popcount is a single machine instruction for obtaining the number of ones in a 32-bit word.
    }
    answer += facdig * fac; // Using base factorial, make next digit be fac
    fac *= (length - i);
    seenletters = seenletters | (1 << digit);
  }
  return answer;
}


// Input: perm, perm's inverse (which needs to be correct in position length - index), perm's length, index, answer = complement of normalization of (index)-prefix of perm, a bitmap named seenpos which should start off at zero for index = 0. 
// Output: bitmap is updated to keep track of the positions in perm of each letter from n - i to n. answer is updated to be complement of normalization of (index + 1)-prefix of perm
void extendnormalizetop(perm_t perm, perm_t inverse, int length, int index, perm_t &answer, uint32_t & seenpos) {
  int i = length - index - 1; 
  int oldpos = getdigit(inverse, i); // position of (length - index) in perm
  int newpos = 0; // will be position of (length - index) in normalization of (i+1)-prefix of perm
  if (oldpos != 0){
    uint32_t temp = seenpos << (32 - oldpos); // Note: shifting by 32 is ill-defined, which is why we explicitly eliminate digit = 0 case.
    newpos = __builtin_popcount(temp);
  }
  answer = setdigit(addpos(answer, newpos), newpos, index);
  seenpos = seenpos | (1 << oldpos);
}

