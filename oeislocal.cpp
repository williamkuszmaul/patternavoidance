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

Oeis ::  Oeis(string filename, int sequencesize, int maxshift)
  : sequencesize(sequencesize)
  , maxshift(maxshift)
{
  
  //cout<<"Building local version of OEIS..."<<endl;
  
  // start by making list of OEIS sequences
  vector < vector < long long > > sequences;
  int maxsequencelength = maxshift + sequencesize - 1; // maximum length of oeis sequence we need to look at
  ifstream input;
  input.open(filename); // file containing list of oeis sequences. Can be downloaded at http://oeis.org/wiki/Welcome#Compressed_Versions
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
  
  // Next build the sequencemap. In particular, for each sequence
  // tempsequence which can be obtained by starting in the <=
  // maxshift-th position of OEIS sequence i,
  // sequencemap[tempsequence] is set to i.
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
      // since earlier oeis sequence tend to be more relevent, better to use them
      if (sequencemap.find(tempsequence) == sequencemap.end()) sequencemap[tempsequence] = i + 1; // to get oeis number
    }
  }
}


 // Given a string containing a sequence separated by spaces of
  // length at least inputshift + sequencesize, extracts the sequence
  // starting with the (inputshift)-th number of line. Indexed so that
  // if inputshift = 0, will start with first entry of line.
Sequence Oeis :: extractusersequence(string line, int inputshift) {
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
  //testsequence.display();
  return testsequence;
}

int Oeis :: getoeisnum(Sequence &sequence) {
  unordered_map<Sequence, int>::const_iterator elt = sequencemap.find(sequence);
  if (elt == sequencemap.end()) return -1;
  return elt->second;
}

// returns false if sequence is detected to start growing as a constant, linear, or quadratic
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

void buildoeisnames(vector <string> & oeisnames) {
  string filename = "names";
  ifstream file(filename);
  string line;
  oeisnames.push_back("No first sequence");
  for (int i = 0; i < 4; i++) getline(file, line);
  while (getline(file, line)) {
    oeisnames.push_back(line);
  }
}

struct oeishitinfo {
  int numhits;
  int oeisnum;
  string additionalinfo;
  oeishitinfo(int a, int b, string c) {
    numhits = a;
    oeisnum = b;
    additionalinfo = c;
  }
  oeishitinfo() {
    numhits = 0;
    oeisnum = -1;
    additionalinfo = "FAILED TO INITIALIZE";
  }
};

// Takes input file and for each line l in input file, writes to output file:
// (1) l (as one line)
// (2) if l does not start with #, checks if l corresponds with an oeis sequence and if so writes the OEIS sequence (as next line)
// If ignoreboring, however, then (2) is NOT printed if the sequence represented by l is detected to be "boring" by allowsequence
void analyzesequencefile(ifstream &inputsequences, ofstream &output, int inputshift, Oeis &OEIS, bool ignoreboring, bool verbose) {
  if (verbose) cout<<"Analyzing sequences..."<<endl;
  vector <string> oeisnames;
  buildoeisnames(oeisnames); 
  


  map<int, oeishitinfo> seensequences;
  string prevline = "";
  string line = "";
  int numwins = 0;
  int numtries = 0;
  int numignored = 0;
  while (getline(inputsequences, line)) {
    output<<line<<endl;
    if (line[0] != '#') { // line is not just commentary
      Sequence testsequence = OEIS.extractusersequence(line, inputshift);
      if (ignoreboring && !allowsequence(testsequence)) {
  	numignored++;
      } else {
	int oeisnum = OEIS.getoeisnum(testsequence);
  	if (oeisnum != -1) {
  	  output<<numtooeis(oeisnum)<<endl;
  	  numwins++;
  	  map<int, oeishitinfo>::const_iterator oe  = seensequences.find(oeisnum);
  	  if (oe == seensequences.end()) {
	    oeishitinfo hitinfo(1, oeisnum, prevline);
  	    seensequences[oeisnum] = hitinfo;
  	  } else {
  	    seensequences[oeisnum].numhits++;
  	  }
  	}
  	numtries++;
      }
    }
    prevline = line;
  }

  // print this regardless of verbose-ness
  cout<<numwins<<" successes out of "<<numtries<<" tries. And "<<numignored<<" ignored."<<endl;
  cout<<"Number distinct sequences: "<<seensequences.size()<<endl;
  
  inputsequences.close();
  output.close();

  if (verbose) {
    std::map<uint64_t, oeishitinfo> tempmap;
    int marker = 0;
    for( std::map<int, oeishitinfo>::iterator iter = seensequences.begin();
	 iter != seensequences.end();
	 ++iter ) {
      cout<<"-----------"<<endl;
      cout<<"Number and number of hits: "<<iter->first<<" "<<iter->second.numhits<<endl; // print out sequences detected in order of oeis number, as well as number of times sequence appears
      cout<<oeisnames[iter->first]<<endl;
      cout<<iter->second.additionalinfo<<endl;
      tempmap[(((uint64_t)(iter->second.numhits))<<22) + marker] = iter->second; // this is a silly hack so that we can in a minute get sequences in order of number of times they appear
      marker++;
    }
    //winfile.close();
    
    cout<<"------------------------------------------------"<<endl;
    
    for( std::map<uint64_t, oeishitinfo>::iterator iter = tempmap.begin();
	 iter != tempmap.end();
	 ++iter ) {
      cout<<"-----------"<<endl;
      cout<<"Number and number of hits: "<<iter->second.oeisnum<<" "<<((iter->first)>>22)<<endl; // print out sequences detected in order of number of times sequence appears
      cout<<oeisnames[iter->second.oeisnum]<<endl;
      cout<<"Example set: "<<iter->second.additionalinfo<<endl;
    }
  }
}

// example usage
// int main(int argc, char* argv[]) {
//   Oeis OEIS("stripped", 8, 15);
//   analyzesequencefile("out-foo", "out-out-foo", 5, OEIS, true);
//   return 0;
// }
