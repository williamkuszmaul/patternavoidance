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


// Given a number, formats it as an OEIS sequence number
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


// good hash function for unsigned long long
//This hash function comes from http://www.concentric.net/~ttwang/tech/inthash.htm
unsigned long long singlehash (unsigned long long key);
//Extension of hash function to sequences of unsigned long longs. Just xors individual hash functions
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
  int maxshift; // consider subsequences starting in positions 1, 2, ..., maxshift
  unordered_map<Sequence, int> sequencemap; // stores pairs (sequence, OEIS number) // Note, only one OEIS number is stored per sequence -- ends up being smallest-valued oeis
  vector <string> oeisnames;
  Oeis(string filename, int sequencesize, int maxshift);
  // Given a string containing a sequence separated by spaces of
  // length at least inputshift + sequencesize, extracts the sequence
  // starting with the (inputshift)-th number of line. Indexed so that
  // if inputshift = 0, will start with first entry of line.
  Sequence extractusersequence(string line, int inputshift);
  // Returns -1 if sequence is not detected in an OEIS sequence. Otherwise, return OEIS number of (arbitrarily selected) OEIS sequence containing this as a subsequence.
  int getoeisnum(Sequence &sequence);
};

// Returns false if sequence is detected to become either constant, quadratic, or cubic. Useful if you want to get rid of boring sequences
bool allowsequence(Sequence &testsequence);

// Takes input file and for each line l in input file, writes to output file:
// (1) l (as one line)
// (2) if l does not start with #, checks if l corresponds with an oeis sequence and if so writes the OEIS sequence (as next line)
// If ignoreboring, however, then (2) is NOT printed if the sequence represented by l is detected to be "boring" by allowsequence
void analyzesequencefile(ifstream &inputsequences, ofstream &output, int inputshift, Oeis &OEIS, bool ignoreboring, bool verbose);

#endif 
