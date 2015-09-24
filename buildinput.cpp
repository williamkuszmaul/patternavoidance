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
#include "buildinput.h"
#include "perm.h"
using namespace std;

// A simple bitmap. Will be used to keep represent a set of patterns
class Bitmap {
public:
  vector <bool> map;
  int size;
  Bitmap(int requestedsize) {
    size = requestedsize;
    map.resize(size);
  }
  void setpos(int pos) {
    map[pos] = true;
  }
  void clearpos(int pos) {
    map[pos] = false;
  }
  bool getpos(int pos) {
    return map[pos];
  }
  void display() {
    for (int i = 0; i < size; i++) cout<<map[i]<<" ";
    cout<<endl;
  }
  int numonbits(){
    int answer = 0;
    for (int i = 0; i < size; i++) answer += (int)map[i];
    return answer;
  }
  void clear() {
    for (int i = 0; i < size; i++) map[i] = false;
  }
};

bool operator> (Bitmap &map1, Bitmap &map2) {
  for (int i = 0; i < map1.size && i < map2.size; i++) {
    if (map1.getpos(i) != map2.getpos(i)) {
      return map1.getpos(i);
    }
  }
  return false;
}

// Permutes a bitmap by sending i-th bit to permutation(i)-th bit
Bitmap shufflebits(Bitmap &bitmap, vector <int> &permutation) {
  Bitmap answer(bitmap.size);
  for (int i = 0; i < bitmap.size; i++) {
    if (bitmap.getpos(i)) answer.setpos(permutation[i]);   
  }
  return answer;
}

// Given a bitmap, returns bitmap where each permutation perm
// represented by a 1 bit is replaced by a
// reverse(inverse(perm)). This corresponds with a rotation if we are
// viewing the trivial-wilf symmetries as the symmetries of a square
Bitmap rotbitmap(Bitmap& bitmap, vector <int> &reverseindices, vector <int> &inverseindices) {
  Bitmap temp = shufflebits(bitmap, inverseindices);
  return shufflebits(temp, reverseindices);
}

// Builds bit maps corresponding with trivially wilf-equivalent pattern sets. Returns whether this bitmap is the lexicographically smallest bitmap of these.
bool isminset(Bitmap& bitmap, vector <int> &reverseindices, vector <int> &inverseindices, vector <int> &complementindices) {
  Bitmap revmap = shufflebits(bitmap, reverseindices); // Corresponds with reflectional symmetry 
  if (bitmap > revmap) return false;
  Bitmap rots1 = bitmap;
  Bitmap rots2 = revmap;
  for (int i = 0; i < 3; i++) {
    rots1 = rotbitmap(rots1, reverseindices, inverseindices);
    rots2 = rotbitmap(rots2, reverseindices, inverseindices);
    if (bitmap > rots1) return false;
    if (bitmap > rots2) return false;
  }
  return true;
}

// Builds all permutations of size finalsize. Adds them to optionvals and the string representation of them to patternoptions
void buildpermutations(perm_t perm, int currentsize, int finalsize, vector <string> &patternoptions, vector <perm_t> &optionvals) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      perm_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize, patternoptions, optionvals);
    }
  }
  if (currentsize == finalsize) {
    string pattern = "";
    for (int i = 0; i < currentsize; i++) pattern = pattern + (char)((char)getdigit(perm, i) + (char)1 + '0');
    patternoptions.push_back(pattern);
    optionvals.push_back(perm);
    //cout<<pattern<<endl;
    //displayperm(perm);
  }
}


// Finds position of val in a vector of values. There is probably a C++ built-in for this...
int find(vector<perm_t> &vec, perm_t val) {
  for (int i = 0; i < vec.size(); i++) {
    if (vec[i] == val) return i;
  }
  assert(1 == 2); // should never get here
  return -1;
}

// Given a bitmap for a set of patterns, checks if map is the unique represetative for its trivial wilf-class.
// If it is, writes the set of patterns to output file.
void entermap(Bitmap & map, ofstream &file, int &numsetstotal, vector <string> &optionsvec, vector <int> &reverseindices, vector <int> &inverseindices, vector <int> &complementindices) { 
  if (isminset(map, reverseindices, inverseindices, complementindices)) {
    vector <string> patternset;
    for (int pos = 0; pos < map.size; pos++) {
      if (map.getpos(pos) == true) patternset.push_back(optionsvec[pos]);
    }
    for (int i = 0; i < patternset.size(); i++) {
      file<<patternset[i];
      if (i < patternset.size() - 1) file<<" ";
    }
    numsetstotal++;
    file<<endl;
  }
}

// replaces map with next bitmap lexicographically to have same number of on bits as this one
// returns false if map is the largest such map
bool getnextmap(Bitmap &map) {
  int i = map.size - 1;
  int numleadingones = 0;
  while (i >= 0 && map.getpos(i) == 1) {
    numleadingones++;
    map.clearpos(i);
    i--;
  }
  while (i >= 0 && map.getpos(i) == 0) {
    i--;
  }
  if (i < 0) return false; 
  map.clearpos(i);
  i++;
  for (int j = 0; j < numleadingones + 1; j++) {
    map.setpos(i+j);
  }
  return true;
}


// Writes to file all pattern sets containing patterns of size
// patternsize. HOWERVER, only writes one pattern set per trivial
// wilf-equivalence class (so that the same computation is not done
// multiple times due to symmetry). This is implemented to be faster
// than repeated calls to the function-interface above.
// Note: only implemented for patternsize < 5. (is not feasible to use for patternset >= 5 anyway)
void writepatternsetstofile(ofstream &file, int patternsize, bool verbose) {
  vector<int> reverseindices; // maintains permutation that permutes pattern options to be the reverse patterns
  vector<int> inverseindices;
  vector<int> complementindices;
  vector <string> options; // stores string representation of S_patternsize
  vector <perm_t> optionvals; // stores each element of S_patternsize as unsigned long long
  buildpermutations(0L, 1, patternsize, options, optionvals); // builds options and optionvals
  int numoptions = options.size();
  for (int i = 0; i < optionvals.size(); i++) { // builds permutations mapping each element of options to its reverse, inverse, and complement respectively
    reverseindices.push_back(find(optionvals, getreverse(optionvals[i], options[i].size())));
    inverseindices.push_back(find(optionvals, getinverse(optionvals[i], options[i].size())));
    complementindices.push_back(find(optionvals, getcomplement(optionvals[i], options[i].size())));
  }

  int numsetstotal = 0;
  Bitmap map(numoptions);
  long long i = 0;
  while (1) {
    if (i >= (1L << numoptions)) break;
    i++;
    map.clear();
    int temp = i;
    int pos = 0;
    while (temp > 0) {
      if (temp%2 == 1) map.setpos(pos); // make bitmap have same bit setup as i
      pos++;
      temp/=2;
    }
    entermap(map, file, numsetstotal, options, reverseindices, inverseindices, complementindices); // output pattern set if it is the unique representative of its trivial wilf-class
  }
  if (verbose) {
    cout<<"Number of sets to be analyzed: "<<numsetstotal<<endl;
  }
}

// Writes to file all pattern sets containing patterns of size
// patternsize. HOWERVER, only writes one pattern set per trivial
// wilf-equivalence class (so that the same computation is not done
// multiple times due to symmetry). This is implemented to be faster
// than repeated calls to the function-interface above.
// Note: only implemented for patternsize < 5. (is not feasible to use for patternset >= 5 anyway)
void writepatternsetstofile(ofstream &file, int setsize, int patternsize, bool verbose) {
  vector<int> reverseindices; // permutation which permutes pattern options to be the reverse patterns
  vector<int> inverseindices;
  vector<int> complementindices;
  vector <string> options; // stores string representation of S_patternsize
  vector <perm_t> optionvals; // stores each element of S_patternsize as unsigned long long
  buildpermutations(0L, 1, patternsize, options, optionvals); // builds options and optionvals
  int numoptions = options.size();
  for (int i = 0; i < optionvals.size(); i++) { // builds permutations mapping each element of options to its reverse, inverse, and complement respectively
    reverseindices.push_back(find(optionvals, getreverse(optionvals[i], options[i].size())));
    inverseindices.push_back(find(optionvals, getinverse(optionvals[i], options[i].size())));
    complementindices.push_back(find(optionvals, getcomplement(optionvals[i], options[i].size())));
  }

  int numsetstotal = 0;
  Bitmap map(numoptions);
  for (int j = 0; j < setsize; j++) map.setpos(j);
  while (1) {
    if(!getnextmap(map)) break;
    entermap(map, file, numsetstotal, options, reverseindices, inverseindices, complementindices);
  }
  if (verbose) {
    cout<<"Number of sets to be analyzed: "<<numsetstotal<<endl;
  }
}

