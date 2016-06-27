#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
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
#include "countavoiders.h"
#include "countpatterns.h"
#include "perm.h"
#include "permutilities.h"
using namespace std;

int main() {
  string patternlist = "123 4231"; // write pattern set as list of strings

  // Example building S_n(Pi)
  vector < list < perm_t > > avoiders; // avoiders[i] will store the avoiders appearing in S_i
  buildavoidersfrompatternlist(patternlist, 8, avoiders); 
  for (int i = 1; i <= 8; i++) {
    cout<<"Number avoiders in S_"<<i<<" of "<<patternlist<<" is "<<avoiders[i].size()<<endl; // print number of avoiders in S_n
    cout<<"Example avoider: "<<endl;
    perm_t perm = *(avoiders[i].begin()); // Permutations are represented in 64-bit integers
    displayperm(perm, i); // Permutations have smallest digit 0 instead of 1 -- sorry for this inconsistancy with the patternset format
    cout<<"In particular, the third digit is "<<getdigit(perm, 2)<<endl; // the function getdigit give you the (i-1)-th digit (i.e., indexing starts at zero for digit positions)
  }

  cout<<"-------------------"<<endl;
  // Example counting |S_n(Pi)| (faster and more memory efficient than actually building)
  vector < uint64_t > numavoiders; // numavoiders[i] will store |S_i(Pi)|
  countavoidersfrompatternlist(patternlist, 10, numavoiders); 
  for (int i = 1; i <= 10; i++) {
    cout<<"Number of permutations in S_"<<i<<" avoiding "<<patternlist<<" is "<<numavoiders[i]<<endl; // print number of avoiders in S_n
  }

  // Example using countpatterns.h
  cout<<"-------------------"<<endl;
  vector < vector <int> > tally;
  vector < vector <int> > completelist;
  countpatterns(patternlist, 11, tally, completelist, false, false);
  // Note that changing the last argument to true would result in
  // completelist not being built, but would also result in a change
  // in the algorithm giving good memory utilization
  cout<<"Using pattern set "<<patternlist<<endl;
  string permtocheck = "52341";
  cout<<"Number of pattern-occurrances in "<<permtocheck<<" is "<<completelist[permtocheck.size()][permtonum(stringtoperm(permtocheck), permtocheck.size())]<<endl;
  permtocheck = "21345";
  cout<<"Number of pattern-occurrances in "<<permtocheck<<" is "<<completelist[permtocheck.size()][permtonum(stringtoperm(permtocheck), permtocheck.size())]<<endl;
  for (int i = 1; i <= 11; i++) {
    cout<<"Number of permutations in S_"<<i<<" with no patterns appearing is "<<tally[i][0]<<endl;
    cout<<"Number of permutations in S_"<<i<<" with 1 pattern appearing is "<<tally[i][1]<<endl;
  }
  return 0;
}
