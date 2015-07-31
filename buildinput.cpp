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

// Note: This seems to be the closest thing to what we're doing:
// http://faculty.valpo.edu/lpudwell/maple/webbook/bookoutput.html

string endfile = "testing.txt";
vector<int> reverseindices; // permutation which permutes pattern options to be the reverse patterns; indexing is same as in 64-bit word. Implemented only for <= 64 pattern options
vector<int> inverseindices;
vector<int> complementindices;


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

inline uint64_t getinverse_local_copy(uint64_t perm, int length) { 
  uint64_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, getdigit(perm, i), i);
  }
  return answer;
}

inline uint64_t getreverse(uint64_t perm, int length) {
  uint64_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, i, getdigit(perm, length - i - 1));
  }
  return answer;
}

inline uint64_t getcomplement(uint64_t perm, int length) {
  uint64_t answer = 0;
  for (int i = 0; i < length; i++) {
    answer = setdigit(answer, i, length - getdigit(perm, i) - 1);
  }
  return answer;
}

// sends i-th bit to permutation(i)-th bit
// uint64_t shufflebits(uint64_t bitmap, vector <int> &permutation) {
//   uint64_t answer = 0;
//   for (int i = 0; i < 64; i++) {
//     if ((bitmap >> i) == 0) return answer; // no more 1-bits left in bit map
//     if ((bitmap & (1L << i)) != 0) answer += (1L << permutation[i]);
//   }
//   return answer;
// }

// sends i-th bit to permutation(i)-th bit
Bitmap shufflebits(Bitmap &bitmap, vector <int> &permutation) {
  Bitmap answer(bitmap.size);
  for (int i = 0; i < bitmap.size; i++) {
    if (bitmap.getpos(i)) answer.setpos(permutation[i]);   
  }
  return answer;
}
  
// uint64_t rotbitmap(uint64_t bitmap) {
//   return shufflebits(shufflebits(bitmap, inverseindices), reverseindices);
// }

Bitmap rotbitmap(Bitmap& bitmap) {
  Bitmap temp = shufflebits(bitmap, inverseindices);
  return shufflebits(temp, reverseindices);
}

bool isminset(Bitmap& bitmap) {
  Bitmap revmap = shufflebits(bitmap, reverseindices);
  if (bitmap > revmap) return false;
  Bitmap rots1 = bitmap;
  Bitmap rots2 = revmap;
  for (int i = 0; i < 3; i++) {
    rots1 = rotbitmap(rots1);
    rots2 = rotbitmap(rots2);
    if (bitmap > rots1) return false;
    if (bitmap > rots2) return false;
  }
  return true;
}

// bool isvalidbitmap(uint64_t map) {
//   //if (__builtin_popcountll(map) == 3 && isminset(map)) return true;
//   if (isminset(map)) return true;
//   return false;
// }

bool isvalidbitmap(Bitmap& map) {
  if (isminset(map)) return true;
  return false;
}

bool isvalidpattern(uint64_t perm, int currentsize) {
  if (currentsize == 5) return true; // check not needed in current implementation
  return false;
}

void buildpermutations(uint64_t perm, int currentsize, int finalsize, vector <string> &patternoptions, vector <uint64_t> &optionvals) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      uint64_t extendedperm = setdigit(addpos(perm, i), i, currentsize);
      buildpermutations(extendedperm, currentsize + 1, finalsize, patternoptions, optionvals);
    }
  }
  if (isvalidpattern(perm, currentsize)) {
    string pattern = "";
    for (int i = 0; i < currentsize; i++) pattern = pattern + (char)((char)getdigit(perm, i) + (char)1 + '0');
    patternoptions.push_back(pattern);
    optionvals.push_back(perm);
    //cout<<pattern<<endl;
    //displayperm(perm);
  }
}


int find(vector<uint64_t> &vec, uint64_t val) {
  for (int i = 0; i < vec.size(); i++) {
    if (vec[i] == val) return i;
  }
  assert(1 == 2);
  return -1;
}

void entermap(Bitmap & map, ofstream &file, int &numsetstotal, vector <string> &optionsvec) { 
  if (isvalidbitmap(map)) {
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


int main(int argc, char* argv[]) {
  
  ofstream file;
  file.open(endfile, std::ofstream::trunc);
  int numoptions = 120;
  int patternsperset = 4;
  int maxpatternsize = 5;
  vector <string> options;
  vector <uint64_t> optionvals;
  buildpermutations(0L, 1, maxpatternsize, options, optionvals);
  
  for (int i = 0; i < optionvals.size(); i++) {
    reverseindices.push_back(find(optionvals, getreverse(optionvals[i], options[i].size())));
    //cout<<reverseindices[i]<<" ";
    inverseindices.push_back(find(optionvals, getinverse_local_copy(optionvals[i], options[i].size())));
    complementindices.push_back(find(optionvals, getcomplement(optionvals[i], options[i].size())));
  }
  //cout<<endl;
  
  int numsetstotal = 0;
  Bitmap map(numoptions);
  for (int j = 0; j < patternsperset; j++) map.setpos(j);
  //int i = 0; // swap out with this commented code to go through all pattern sets of each size
  while (1) {
    if(!getnextmap(map)) break;
    // if (i >= (1L << numoptions)) break;
    // i++;
    // map.clear();
    // int temp = i;
    // int pos = 0;
    // while (temp > 0) {
    //   if (temp%2 == 1) map.setpos(pos);
    //   pos++;
    //   temp/=2;
    // }
    entermap(map, file, numsetstotal, options);
  }
  cout<<"Number of sets to be analyzed: "<<numsetstotal<<endl;
  cout<<"Output is in file: "<<endfile<<endl;
  file.close();
  return 0;
}
