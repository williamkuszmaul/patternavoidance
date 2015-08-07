#ifndef _OEISLOCAL_H 
#define _OEISLOCAL_H // To avoid header being included twice in complilation process.


#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <algorithm>
#include <time.h> 
#include <stdlib.h>
#include <bitset>
#include <string>
#include <vector>
#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <queue>
#include "hashdb.h"
#include <sys/time.h>
#include "fastavoidance.h"
using namespace std;


// helper functions ---------------------------------------------------
string numtooeis(int num);

// Sequence class -------------------------------------------------
class Sequence {
public:
  vector <long long> data;
  Sequence(int size);
  Sequence(vector <long long> vec);
  bool operator==(const Sequence &other) const
  {
    return (data == other.data);
  }
  void display();
};


unsigned long long singlehash (unsigned long long key);

namespace std {
  template <>
  struct hash<Sequence>
  {
    std::size_t operator()(const Sequence & key) const
    {
      uint64_t answer = 0;
      for (int i = 0; i < key.data.size(); i++) answer = answer ^ singlehash(key.data[i]);
      return answer;
    }
  };
}

// OEIS class ---------------------------------------------------------------------

class Oeis {
public:
  int sequencesize; // sequence size
  int maxshift; // consider subsequences starting in pos <= maxshift-th
  unordered_map<Sequence, int> sequencemap; // stores pairs (sequence, OEIS number) // Note, only one OEIS number is stored per sequence -- ends up being smallest-valued oeis
  Oeis(string filename, int sequencesize, int maxshift);
  Sequence extractusersequence(string line, int inputshift);
  int getoeisnum(Sequence &sequence);
};

bool allowsequence(Sequence &testsequence);
void analyzesequencefile(ifstream &inputsequences, ofstream &output, int inputshift, Oeis &OEIS, bool verbose);

#endif 
