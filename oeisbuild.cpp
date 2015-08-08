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
using namespace std;

const int sequencesize = 9; // sequence size
const int maxshift = 15; // consider subsequences starting in pos <= maxshift-th
const int inputshift= 3;
bool ignoreon = true;

template <int SEQUENCESIZE> 
class Sequence {
  public:
    unsigned int data[SEQUENCESIZE];
    Sequence() {}
    Sequence(int *array) {
      for (int i = 0 ; i < SEQUENCESIZE; i++) {
	data[i] = array[i];
    }
  }
  bool operator==(const Sequence &other) const
  {
    for (int i = 0; i < SEQUENCESIZE; i++) {
      if (data[i] != other.data[i]) return false;
    }
    return true;
  }
};


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


namespace std {
template <>
struct hash<Sequence<sequencesize>>
{
  std::size_t operator()(const Sequence<sequencesize> & key) const
  {
    uint64_t answer = 0;
    for (int i = 0; i < sequencesize; i++) answer = answer ^ singlehash(key.data[i]);
    return answer;
  }
};
}

string numtooeis(int num) {
  string answer = "";
  for (int i = 0; i < 6; i++) {
    answer = (char)((char)(num % 10) + '0')+ answer;
    num /= 10;
  }
  return "A" + answer;
}

bool allowsequence(Sequence<sequencesize> &testsequence) {
  Sequence<sequencesize> deriv1;
  Sequence<sequencesize> deriv2;
  Sequence<sequencesize> deriv3;
  if (ignoreon) {
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
  return true;
}

int main(int argc, char* argv[]) {
  cout<<sequencesize<<" is length of analysis, starting in "<<(inputshift + 1)<<"-th position of input sequences, meaning inputsequences must be length at least "<<sequencesize+inputshift<<" and maxshift is "<<maxshift<<endl; 
  int maxsequencelength = maxshift + sequencesize - 1; // maximum length of oeis sequence we need to look at
  vector < vector <int> > sequences;

  if (argc != 4) {
    cout<<"Program takes argument of oeis file name (in stripped form obtainable at http://oeis.org/wiki/Welcome#Compressed_Versions  and file containing sequences of interest and file for output of sequences which worked"<<endl;
    return 0;
  }  
  cout<<"Building local version of OEIS..."<<endl;
  string filename = argv[1];
  ifstream input;
  input.open(filename);
  string line;
  for (int i = 0; i < 4; i++) getline(input, line); // get rid of garbage
  while (getline(input, line)) {
    vector <int> sequence;
    int index = 9;
    int nextval = 0;
    while (sequence.size() < maxsequencelength && index < line.size() - 1) {
      if (line[index] == ',') {
	sequence.push_back(nextval);
	//cout<<nextval<<" ";
	nextval = 0;
      } else {
	nextval = nextval * 10 + (int)(line[index] - '0');
      }
      index++;
    }
    // cout<<endl;
    sequences.push_back(sequence);
    //cout<<line<<endl;
  }
  input.close();


  unordered_map<Sequence<sequencesize>, int> sequencemap; // stores pairs (sequence, OEIS number) // Note, only one OEIS number is stored per sequence -- ends up being smallest-valued oe
  for (int i = 0; i < sequences.size(); i++) {
    for (int start = 0; start < maxshift; start++) { // start pos
      Sequence<sequencesize> tempsequence;
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

  // Example usage:
   // int testdata[] = {4, 8, 13, 25, 44, 83, 152, 286, 538, 1020};
   // int testdata2[] = {4, 8, 13, 25, 44, 83, 152, 286, 11, 1020};
   // Sequence<sequencesize> testsequence(testdata);
   // Sequence<sequencesize> testsequence2(testdata2);
   // cout<<(sequencemap.find(testdata) == sequencemap.end())<<endl;
   // cout<<(sequencemap.find(testdata2) == sequencemap.end())<<endl;

  cout<<"Analyzing sequences..."<<endl;
  string filename2 = argv[2];
  string winfilename = argv[3];
  ofstream winfile;
  winfile.open(winfilename, std::ofstream::trunc);


  ifstream inputsequences;
  inputsequences.open(filename2);
  string outputfile = "out-"+filename2;
  ofstream output;
  output.open(outputfile, std::ofstream::trunc);

  map<int, int> seensequences;

  string prevline = "";
  line = "";
  int numwins = 0;
  int numtries = 0;
  int numignored = 0;
  while (getline(inputsequences, line)) {
    output<<line<<endl;
    if (line[0] != '#') { // line is not just commentary
      Sequence<sequencesize> testsequence;
      int index = 0;
      int nextval = 0;
      int sequencepos = 0;
      while (sequencepos < inputshift + sequencesize) {
	assert(index < line.size());
	if (line[index] == ' ') {
	  if (sequencepos >= inputshift) testsequence.data[sequencepos - inputshift] = nextval;
	  //cout<<nextval<<" ";
	  nextval = 0;
	  sequencepos++;
	} else {
	  nextval = nextval * 10 + (int)(line[index] - '0');
	}
	index++;
      }
      //cout<<endl;
      if (!allowsequence(testsequence)) {
	numignored++;
      } else {
	unordered_map<Sequence<sequencesize>, int>::const_iterator elt = sequencemap.find(testsequence);
	if (elt != sequencemap.end()) {
	  // winfile<<prevline.substr(1)<<endl; // get rid of hashtag at start
	  output<<numtooeis(elt->second)<<endl;
	  numwins++;
	  //if (seensequences.find(elt->second) == seensequences.end()) cout<<numtooeis(elt->second)<<endl;
	  map<int, int>::const_iterator oe  = seensequences.find(elt->second);
	  if (oe == seensequences.end()) {
	    seensequences[elt->second] = 1;
	  } else {
	    seensequences[elt->second]++;
	  }
	}
	numtries++;
      }
    }
    prevline = line;
  }
  cout<<numwins<<" successes out of "<<numtries<<" tries. And "<<numignored<<" ignored."<<endl;
  cout<<"Number distinct sequences: "<<seensequences.size()<<endl;
  cout<<"Output is in file: "<<outputfile<<endl;
  inputsequences.close();
  output.close();

  std::map<uint64_t, int> tempmap;
  int marker = 0;
  for( std::map<int, int>::iterator iter = seensequences.begin();
     iter != seensequences.end();
     ++iter ) {
    winfile<<iter->first<<" ";
    cout<<iter->first<<" "<<iter->second<<endl;
    tempmap[(((uint64_t)(iter->second))<<22) + marker] = iter->first;
    marker++;
  }
  winfile.close();

  cout<<"------------------------------------------------"<<endl;

  for( std::map<uint64_t, int>::iterator iter = tempmap.begin();
       iter != tempmap.end();
       ++iter ) {
    cout<<iter->second<<" "<<((iter->first)>>22)<<endl;
  }

  return 0;
}
