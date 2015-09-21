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
#include <list>
// #include <unordered_set>
#include <queue>
#include "hashdb.h"
#include <sys/time.h>
#include "fastavoidance.h"
#include "countpatterns.h"
#include "perm.h"
using namespace std;

int main() {
  string patternlist = "123 213"; // write pattern set as list of strings

  // example building S_n(Pi)
  vector < list < uint64_t > > avoiders; // avoiders[i] will store the avoiders appearing in S_i
  buildavoidersfrompatternlist(patternlist, 8, avoiders); 
  for (int i = 1; i <= 8; i++) {
    cout<<"Number avoiders in S_"<<i<<" of "<<patternlist<<" is "<<avoiders[i].size()<<endl; // print number of avoiders in S_n
    cout<<"Example avoider: "<<endl;
    uint64_t perm = *(avoiders[i].begin()); // Permutations are represented in 64-bit integers
    displayperm(perm, i); // Permutations have smallest digit 0 instead of 1 -- sorry for this inconsistancy with the patternset format
    cout<<"In particular, the third digit is "<<getdigit(perm, 2)<<endl; // the function getdigit give you the (i-1)-th digit (i.e., indexing starts at zero for digit positions)
  }

  cout<<"-------------------"<<endl;
  // example counting |S_n(Pi)| (slightly faster than building)
  patternlist = "123 213"; // write pattern set as list of strings
  vector < uint64_t > numavoiders; // numavoiders[i] will store |S_i(Pi)|
  countavoidersfrompatternlist(patternlist, 10, numavoiders); 
  for (int i = 1; i <= 10; i++) {
    cout<<"Number of permutations in S_"<<i<<" avoiding "<<patternlist<<" is "<<numavoiders[i]<<endl; // print number of avoiders in S_n
  }

  // example using countpatterns
  cout<<"-------------------"<<endl;
  patternlist = "123 213";
  vector < vector <int> > tally;
  vector < vector <int> > completelist;
  countpatterns(patternlist, 10, tally, completelist, false, false);
  string permtocheck = "21345";
  cout<<"Using pattern set "<<patternlist<<endl;
  cout<<"Number of pattern-occurrances in "<<permtocheck<<" is "<<completelist[permtocheck.size()][permtonum(stringtoperm(permtocheck), permtocheck.size())]<<endl;
  cout<<"Number of permutations in S_10 with no patterns appearing is "<<tally[10][0]<<endl;
  cout<<"Number of permutations in S_10 with 1 pattern appearing is "<<tally[10][1]<<endl;
  return 0;
}
