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
using namespace std;

string endfile = "input.txt";

bool isvalidbitmap(uint64_t map) {
  if (__builtin_popcountll(map) == 4) return true;
  return false;
}

bool isvalidpattern(uint64_t perm, int currentsize) {
  if (currentsize == 4) return true;
  return false;
}

void buildpermutations(uint64_t perm, int currentsize, int finalsize, vector <string> &patternoptions) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize, patternoptions);
    }
  }
  if (isvalidpattern(perm, currentsize)) {
    string pattern = "";
    for (int i = 0; i < currentsize; i++) pattern = pattern + (char)((char)getdigit(perm, i) + (char)1 + '0');
    patternoptions.push_back(pattern);
  }
}


int main(int argc, char* argv[]) {
  
  ofstream file;
  file.open(endfile, std::ofstream::trunc);
  int numoptions = 24;
  vector <string> options;
  buildpermutations(0L, 1, 4, options);
  for (uint64_t i = 0; i < (1L << numoptions); i++) {
    if (isvalidbitmap(i)) {
      vector <string> patternset;
      uint64_t bitmap = i;
      int index = 0;
      while (bitmap > 0) {
	if (bitmap%2 == 1) patternset.push_back(options[index]);
	bitmap /= 2;
	index++;
      }
      for (int i = 0; i < patternset.size(); i++) {
	file<<patternset[i];
	if (i < patternset.size() - 1) file<<" ";
      }
      file<<endl;
    }
  }
  file.close();
  return 0;
}
