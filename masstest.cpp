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
#include "oeislocal.h"
using namespace std;

void setstosequences(string inputfile, string outputfile, int maxpermsize) {
  ifstream setsfilein;
  setsfilein.open(inputfile);
  ofstream sequencesfile;
  sequencesfile.open(outputfile, std::ofstream::trunc);
  countavoidersfromfile(setsfilein, sequencesfile, maxpermsize, false);
  setsfilein.close();
  sequencesfile.close();
}

void buildwilfeqsets(string setsfilename, int patternsize) {
  cout<<"Building sets of "<<patternsize<<"-patterns up to trivial Wilf equivalence..."<<endl;
  ofstream setsfile;
  setsfile.open(setsfilename, std::ofstream::trunc);
  writepatternsetstofile(setsfile, patternsize, false);
  setsfile.close();
}

void analyzefirstsequenceset(Oeis &OEIS, string sequencesfilename, string setsfilename2, int minpermsize, int minrevisedsetsize, unordered_map<int, string> &oeismatchesbynum) {
  int firstzeropos = 4, secondzeropos = 3, thirdzeropos = 4;
  ifstream sequencesfilein;
  sequencesfilein.open(sequencesfilename);
  ofstream nextlevelsets;
  nextlevelsets.open(setsfilename2, std::ofstream::trunc);
  vector <patternsetinfo> matches;
  int numattempts = 0;
  fillpatternsetinfo(sequencesfilein, OEIS, minpermsize - 1, matches, numattempts); // start with the minpermsize-th entry
  sequencesfilein.close();
  int constmatches = 0;
  int linmatches = 0;
  int quadmatches = 0;
  unordered_set<int> consthits;
  unordered_set<int> linhits;
  unordered_set<int> quadhits;
  for (int i = 0; i < matches.size(); i++) {
    patternsetinfo set = matches[i];
    if (iszeroby(ithderivative(set.sequence, 1), firstzeropos)) { // look at first derivative starting in 4-th position of sequence (since some sequence become constant late
      consthits.insert(set.oeisnum);
      constmatches++;
    } else {
      if (iszeroby(ithderivative(set.sequence, 2), secondzeropos)) { // look at second derivative starting with 3rd term
	linmatches++;
	linhits.insert(set.oeisnum);
      } else {
	if (iszeroby(ithderivative(set.sequence, 3), thirdzeropos)) { // look at third derivative starting with 4th term
	  quadmatches ++;
	  quadhits.insert(set.oeisnum);
	} else {
	  //	  cout<<"entry "<<set.oeisnum<<" "<<set.patternset<<endl;
	  if (oeismatchesbynum.find(set.oeisnum) == oeismatchesbynum.end()) {
	    oeismatchesbynum[set.oeisnum] = set.patternset;
	  }
	  else {
	    oeismatchesbynum[set.oeisnum] = oeismatchesbynum[set.oeisnum] + "\n" + set.patternset;
	  }
	  if (std::count(set.patternset.begin(), set.patternset.end(), ' ') >= minrevisedsetsize) nextlevelsets<<set.patternset<<endl;
	}
      }
    }
  }
  cout<<"Total of "<<matches.size()<<" matches out of "<<numattempts<<" attempts"<<endl;
  cout<<constmatches<<" matches constant by first derivative position "<<firstzeropos<<endl;
  cout<<linmatches<<" of remaining matches linear by second derivative position "<<secondzeropos<<endl;
  cout<<quadmatches<<" of remaining matches quadratic by third derivative  position "<<secondzeropos<<endl;
  cout<<"Number distinct distinct oeis sequences matching to constant sequences: "<<consthits.size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over linear sequences: "<<linhits.size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over quadratic sequences: "<<quadhits.size()<<endl;
  cout<<"Number distinct oeis sequences for remaining: "<<oeismatchesbynum.size()<<endl;
  // for (std::unordered_map<int, string>::iterator iter = oeismatchesbynum.begin();
  //      iter != oeismatchesbynum.end();
  //      ++iter) {
  //   //cout<<"second type of entry"<<iter->first<<" "<<iter->second<<endl;
  // }
  nextlevelsets.close();
}

void analyzesecondsequenceset(Oeis &OEIS, string sequencesfilename, int minpermsize, unordered_map<int, string> &oeismatchesbynum) {
  ifstream sequencesfilein;
  sequencesfilein.open(sequencesfilename);
  vector <patternsetinfo> matches;
  int numattempts = 0;
  fillpatternsetinfo(sequencesfilein, OEIS, minpermsize - 1, matches, numattempts); // start with the minpermsize-th entry
  sequencesfilein.close();
  for (int i = 0; i < matches.size(); i++) {
    patternsetinfo set = matches[i];
    if (oeismatchesbynum.find(set.oeisnum) == oeismatchesbynum.end()) {
      oeismatchesbynum[set.oeisnum] = set.patternset;
    }
    else {
      oeismatchesbynum[set.oeisnum] = oeismatchesbynum[set.oeisnum] + "\n" + set.patternset;
    }
  }
  cout<<"Total of "<<matches.size()<<" oeis matches for filtered sequences"<<endl;
  cout<<"Total of "<<oeismatchesbynum.size()<<" distinct oeis matches for filtered sequences"<<endl;
}

int main() {
  // Writes all pattern sets of S_4 up to trivial wilf-equivalence to setsfilename
  // For each one, computes |S_1(Pi)|, ..., |S_maxpermsize(Pi)| and writes to sequencesfilename
  // For each sequence in sequencesfilename, looks up S_minpermsize(Pi), ..., S_maxpermsize(Pi) in OEIS and writes found oeis sequences to oeismatchesfilename
  
  int maxpermsize = 13; //
  int minpermsize = 5;
  int patternsize = 3;
  int revisedmaxpermsize = 15;
  int minrevisedsetsize = 2;
  string setsfilename = "testallfours-sets.txt";
  string sequencesfilename = "testallfours-out13";
  string setsfilename2 = "testallfours-revisedsets";
  string sequencesfilename2 = "testallfours-revisedout15.txt";

  buildwilfeqsets(setsfilename, patternsize);
  cout<<"Building sequences up to "<<maxpermsize<<endl;
  setstosequences(setsfilename, sequencesfilename, maxpermsize);
  
  cout<<"Building local version of OEIS..."<<endl;
  Oeis OEIS("/data/williamkuszmaul/stripped", "/data/williamkuszmaul/names", maxpermsize - minpermsize + 1, 15); // Note: we allow sequences to start in any of positions 1, 2, ..., 15
  cout<<"Continuing analysis"<<endl;

  unordered_map<int, string> oeismatchesbynum;
  analyzefirstsequenceset(OEIS, sequencesfilename, setsfilename2, minpermsize, minrevisedsetsize, oeismatchesbynum);
  cout<<"Building sequences for filtered sets up to "<<revisedmaxpermsize<<endl;
  setstosequences(setsfilename2, sequencesfilename2, revisedmaxpermsize);
  unordered_map<int, string> oeismatchesbynum2;
  cout<<"Analyzing filtered sequences"<<endl;
  analyzesecondsequenceset(OEIS, sequencesfilename2, minpermsize, oeismatchesbynum2);
  
  return 0;  
}
