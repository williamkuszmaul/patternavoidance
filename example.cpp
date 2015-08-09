// To compile, use: g++ -O3 -std=c++11 -march=native -g -o example fastavoidance.cpp hashdb.cpp example.cpp
// To run, use: ./example

#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
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
#include "perm.h"
using namespace std;

int main() {
  string patternset = "2431 4231 1432 4132"; // write pattern set as list of strings

  // example building S_n(Pi)
  vector < vector < uint64_t > > avoiders; // avoiders[i] will store the avoiders appearing in S_i
  buildavoidersfrompatternlist(patternset, 8, avoiders); 
  for (int i = 1; i <= 8; i++) {
    cout<<avoiders[i].size()<<endl; // print number of avoiders in S_n
    cout<<"Example avoider: "<<endl;
    uint64_t perm = avoiders[i][0]; // Permutations are represented in 64-bit integers
    displayperm(perm, i); // Permutations have smallest digit 0 instead of 1 -- sorry for this inconsistancy with the patternset format
    cout<<"In particular, the third digit is "<<getdigit(perm, 2)<<endl; // the function getdigit give you the (i-1)-th digit (i.e., indexing starts at zero for digit positions)
  }

  // example counting |S_n(Pi)| (slightly faster than building)
  vector < uint64_t > numavoiders; // numavoiders[i] will store |S_i(Pi)|
  countavoidersfrompatternlist(patternset, 8, numavoiders); 
  for (int i = 1; i <= 8; i++) {
    cout<<numavoiders[i]<<endl; // print number of avoiders in S_n
  }

  return 0;
}
