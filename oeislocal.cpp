// First argument is file in format of stripped file found at http://oeis.org/wiki/Welcome#Compressed_Versions
// Second argument is file containing one sequence per line, each of length sequencesize (const global variable -- that way I can hardcode structs to have it as their size).
// Sequences in second file are separated by spaces. There is no space at start, but there NEEDS TO BE a leading space at the end

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
#include "oeislocal.h"
using namespace std;

// helper functions ---------------------------------------------------
string numtooeis(int num) {
  string answer = "";
  for (int i = 0; i < 6; i++) {
    answer = (char)((char)(num % 10) + '0')+ answer;
    num /= 10;
  }
  return "A" + answer;
}

// Sequence class -------------------------------------------------

Sequence :: Sequence(int size) {
  data.resize(size);
}

Sequence :: Sequence(vector <long long> vec) {
  data = vec;
}

void Sequence :: display() {
  for (int i = 0; i < data.size(); i++) {
    cout<<data[i]<<" ";
  }
  cout<<endl;
}

unsigned long long singlehash (unsigned long long key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key; //assume maxsize is power of two
  //This hash function comes from http://www.concentric.net/~ttwang/tech/inthash.htm
}


// OEIS class ---------------------------------------------------------------------

Oeis ::  Oeis(string filename, int sequencesize, int maxshift, int inputshift)
  : sequencesize(sequencesize)
  , maxshift(maxshift)
  , inputshift(inputshift)
{
  
  cout<<"Building local version of OEIS..."<<endl;
  
  // start by making list of OEIS sequences
  vector < vector < long long > > sequences;
  int maxsequencelength = maxshift + sequencesize - 1; // maximum length of oeis sequence we need to look at
  ifstream input;
  input.open(filename);
  string line;
  for (int i = 0; i < 4; i++) getline(input, line); // get rid of garbage
  while (getline(input, line)) {
    vector < long long > sequence;
    int index = 9;
    long long nextval = 0;
    while (sequence.size() < maxsequencelength && index < line.size() - 1) {
      if (line[index] == ',') {
	sequence.push_back(nextval);
	nextval = 0;
      } else {
	nextval = nextval * 10 + (int)(line[index] - '0');
      }
      index++;
    }
    sequences.push_back(sequence);
  }
  input.close();
  
  // Next build the hashtables
  for (int i = 0; i < sequences.size(); i++) {
    for (int start = 0; start < maxshift; start++) { // start pos
      Sequence tempsequence(sequencesize);
      for (int place = 0; place < sequencesize; place++) {
	if (start + place < sequences[i].size()) {
	  tempsequence.data[place] = sequences[i][start + place];
	} else {
	  tempsequence.data[place] = 0; // if we go off the end of a sequence, just fill in rest with zeros
	}
      }
      sequencemap[tempsequence] = i + 1; // to get oeis number
    }
  }
}

Sequence Oeis :: extractusersequence(string line) {
  Sequence testsequence(sequencesize);
  int index = 0;
  int nextval = 0;
  int sequencepos = 0;
  while (sequencepos < inputshift + sequencesize) {
    assert(index < line.size());
    if (line[index] == ' ') {
      if (sequencepos >= inputshift) testsequence.data[sequencepos - inputshift] = nextval;
      nextval = 0;
      sequencepos++;
      } else {
      nextval = nextval * 10 + (int)(line[index] - '0');
    }
    index++;
  }
  testsequence.display();
  return testsequence;
}

int Oeis :: getoeisnum(Sequence &sequence) {
  unordered_map<Sequence, int>::const_iterator elt = sequencemap.find(sequence);
  if (elt == sequencemap.end()) return -1;
  return elt->second;
}

bool allowsequence(Sequence &testsequence) {
  uint64_t sequencesize = testsequence.data.size();
  Sequence deriv1(sequencesize);
  Sequence deriv2(sequencesize);
  Sequence deriv3(sequencesize);
  for (int i = 1; i < sequencesize; i++) {
    deriv1.data[i] = testsequence.data[i] - testsequence.data[i - 1];
  }
  for (int i = 2; i < sequencesize; i++) {
    deriv2.data[i] = deriv1.data[i] - deriv1.data[i - 1];
  }
  for (int i = 3; i < sequencesize; i++) {
    deriv3.data[i] = deriv2.data[i] - deriv2.data[i - 1];
  }
  
  bool cont = false;
  for (int i = 4; i < sequencesize - 1; i++) { // could start i as low as 1, which would allow more sequences to slip through
    if (deriv1.data[i] != deriv1.data[i+1]) cont = true;
  }
  if (!cont) return false;
  cont = false;
  for (int i = 3; i < sequencesize - 1; i++) { // could start i as low as 2
    if (deriv2.data[i] != deriv2.data[i+1]) cont = true;
  }
  if (!cont) return false;
  cont = false;
  for (int i = 4; i < sequencesize - 1; i++) { // could start i as low as 3
    if (deriv3.data[i] != deriv3.data[i+1]) cont = true;
  }
  if (!cont) return false;
  return true;
}

void analyzesequencefile(string infile, string outfile, Oeis &OEIS, bool verbose) {
  if (verbose) cout<<"Analyzing sequences..."<<endl;
  // string winfilename = argv[3];
  // ofstream winfile;
  // winfile.open(winfilename, std::ofstream::trunc);

  ifstream inputsequences;
  inputsequences.open(infile);
  ofstream output;
  output.open(outfile, std::ofstream::trunc);

  map<int, int> seensequences;
  string prevline = "";
  string line = "";
  int numwins = 0;
  int numtries = 0;
  int numignored = 0;
  while (getline(inputsequences, line)) {
    output<<line<<endl;
    if (line[0] != '#') { // line is not just commentary
      Sequence testsequence = OEIS.extractusersequence(line);
      if (!allowsequence(testsequence)) {
  	numignored++;
      } else {
	int oeisnum = OEIS.getoeisnum(testsequence);
  	if (oeisnum != -1) {
  	  output<<numtooeis(oeisnum)<<endl;
  	  numwins++;
  	  map<int, int>::const_iterator oe  = seensequences.find(oeisnum);
  	  if (oe == seensequences.end()) {
  	    seensequences[oeisnum] = 1;
  	  } else {
  	    seensequences[oeisnum]++;
  	  }
  	}
  	numtries++;
      }
    }
    prevline = line;
  }
  if (verbose) {
    cout<<numwins<<" successes out of "<<numtries<<" tries. And "<<numignored<<" ignored."<<endl;
    cout<<"Number distinct sequences: "<<seensequences.size()<<endl;
    cout<<"Output is in file: "<<outfile<<endl;
  }
  inputsequences.close();
  output.close();

  if (verbose) {
    std::map<uint64_t, int> tempmap;
    int marker = 0;
    for( std::map<int, int>::iterator iter = seensequences.begin();
	 iter != seensequences.end();
	 ++iter ) {
      //  winfile<<iter->first<<" ";
      cout<<iter->first<<" "<<iter->second<<endl;
      tempmap[(((uint64_t)(iter->second))<<22) + marker] = iter->first;
      marker++;
    }
    //winfile.close();
    
    cout<<"------------------------------------------------"<<endl;
    
    for( std::map<uint64_t, int>::iterator iter = tempmap.begin();
	 iter != tempmap.end();
	 ++iter ) {
      cout<<iter->second<<" "<<((iter->first)>>22)<<endl;
    }
  }
}

int main(int argc, char* argv[]) {
  Oeis OEIS("stripped", 8, 15, 5);
  analyzesequencefile("out-foo", "out-out-foo", OEIS, true);
  return 0;
}
