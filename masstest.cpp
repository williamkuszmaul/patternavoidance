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
  countavoidersfromfile_parallel(setsfilein, sequencesfile, maxpermsize);
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
  ifstream sequencesfilein;
  sequencesfilein.open(sequencesfilename);
  ofstream nextlevelsets;
  nextlevelsets.open(setsfilename2, std::ofstream::trunc);
  vector <patternsetinfo> matches;
  int numattempts = 0, numdistinctattempts = 0;
  fillpatternsetinfo(sequencesfilein, OEIS, minpermsize - 1, matches, numattempts, numdistinctattempts); // start with the minpermsize-th entry
  sequencesfilein.close();
  vector <int> numhitswithdeg(4);
  unordered_set<int> hitswithdeg[4];
  for (int i = 0; i < matches.size(); i++) {
    patternsetinfo set = matches[i];
    bool cont = true;
    int degree = 4;
    if (iszeroby(ithderivative(set.sequence, 4), 5)) degree = 3;
    if (iszeroby(ithderivative(set.sequence, 3), 5)) degree = 2;
    if (iszeroby(ithderivative(set.sequence, 2), 5)) degree = 1;
    if (iszeroby(ithderivative(set.sequence, 1), 5)) degree = 0;
    if (degree < 4) {
      numhitswithdeg[degree]++;
      hitswithdeg[degree].insert(set.oeisnum);
    } else {
      if (oeismatchesbynum.find(set.oeisnum) == oeismatchesbynum.end()) {
	oeismatchesbynum[set.oeisnum] = set.patternset;
      } else {
	oeismatchesbynum[set.oeisnum] = oeismatchesbynum[set.oeisnum] + "\n" + set.patternset;
      }
      if (std::count(set.patternset.begin(), set.patternset.end(), ' ') >= minrevisedsetsize)  nextlevelsets<<set.patternset<<endl;
    }
  }
  cout<<"Total of "<<matches.size()<<" matches out of "<<numattempts<<" attempts, with what appears to be "<<numdistinctattempts<<" total Wilf-classes"<<endl;
  cout<<numhitswithdeg[0]<<" matches constant by first derivative position 5"<<endl;
  cout<<numhitswithdeg[1]<<" of remaining matches linear by second derivative position 5"<<endl;
  cout<<numhitswithdeg[2]<<" of remaining matches quadratic by second derivative position 5"<<endl;
  cout<<numhitswithdeg[3]<<" of remaining matches cubic by second derivative position 5"<<endl;
  cout<<"Number distinct distinct oeis sequences matching to constant sequences: "<<hitswithdeg[0].size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over linear sequences: "<<hitswithdeg[1].size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over quadratic sequences: "<<hitswithdeg[2].size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over cubic sequences: "<<hitswithdeg[3].size()<<endl;
  cout<<"Number distinct oeis sequences for remaining: "<<oeismatchesbynum.size()<<endl;
  nextlevelsets.close();
}

void analyzesecondsequenceset(Oeis &OEIS, string sequencesfilename, int minpermsize, unordered_map<int, string> &oeismatchesbynum) {
  ifstream sequencesfilein;
  sequencesfilein.open(sequencesfilename);
  vector <patternsetinfo> matches;
  int numattempts = 0, numdistinctattempts = 0;
  fillpatternsetinfo(sequencesfilein, OEIS, minpermsize - 1, matches, numattempts, numdistinctattempts); // start with the minpermsize-th entry
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

void  compareoeismatches( unordered_map<int, string> &oeismatchesbynum,  unordered_map<int, string> &oeismatchesbynum2) {
  cout<<"Elements appearing only in first analysis of revised data set: "<<endl;
  for (std::unordered_map<int, string>::iterator iter = oeismatchesbynum.begin();
       iter != oeismatchesbynum.end();
       ++iter) {
    int numhits = std::count(iter->second.begin(), iter->second.end(), '\n') + 1;
    if (oeismatchesbynum2.find(iter->first) == oeismatchesbynum2.end() || oeismatchesbynum2[iter->first] != oeismatchesbynum[iter->first]) {
      cout<<numtooeis(iter->first)<<endl<<numhits<<endl;
      //cout<<iter->second<<endl;
    }
  }
  cout<<"Elements appearing only in second analysis of revised data set: "<<endl;
  for (std::unordered_map<int, string>::iterator iter = oeismatchesbynum2.begin();
       iter != oeismatchesbynum2.end();
       ++iter) {
    int numhits = std::count(iter->second.begin(), iter->second.end(), '\n') + 1;
    if (oeismatchesbynum.find(iter->first) == oeismatchesbynum.end() || oeismatchesbynum[iter->first] != oeismatchesbynum2[iter->first]) {
      cout<<numtooeis(iter->first)<<endl<<numhits<<endl;
      //cout<<iter->second<<endl;
    }
  }
}

void  displaymatches(Oeis &OEIS, unordered_map<int, string> &oeismatchesbynum3, string revisedall, bool verbose) {
  ofstream output;
  output.open(revisedall, std::ofstream::trunc);
   for (std::unordered_map<int, string>::iterator iter = oeismatchesbynum3.begin();
       iter != oeismatchesbynum3.end();
       ++iter) {
     int numhits = std::count(iter->second.begin(), iter->second.end(), '\n') + 1;
     output<<OEIS.oeisnames[iter->first]<<endl;
     output<<numhits<<" times"<<endl;
     if (verbose) output<<iter->second<<endl;
   }

}

int main() {
  bool data_built = true; // To just run OEIS analysis since sequence files prebuilt
  int maxpermsize = 13; 
  int minpermsize = 5;
  int patternsize = 4;
  int revisedmaxpermsize = 15;
  int minrevisedsetsize = 5;
  string setsfilename = "testallfours-sets.txt";
  string sequencesfilename = "testallfours-out13TEMP";
  string setsfilename2 = "testallfours-revisedsets";
  string sequencesfilename2 = "testallfours-revisedout15.txt";
  string revisedall = "testallfours-revised15final.txt";
  string revisedbrief = "testallfours-revised15final-brief.txt";

  if (!data_built) buildwilfeqsets(setsfilename, patternsize);
  cout<<"Building sequences up to "<<maxpermsize<<"..."<<endl;
  if (!data_built) setstosequences(setsfilename, sequencesfilename, maxpermsize); // this computation takes around 20 minutes with our code, about 2 hours 20 minutes with permlab. Haven't had chance to run the second computation with permlab. Our algorithm would have do even better (relative to permlab) for larger sets of larger patterns, or larger n.
  cout<<"Building local version of OEIS..."<<endl;
  Oeis OEIS("/data/williamkuszmaul/stripped", "/data/williamkuszmaul/names", maxpermsize - minpermsize + 1, 15); //Note: we allow sequences to start in any of positions 1, 2, ..., 15
  cout<<"Continuing analysis..."<<endl;

  unordered_map<int, string> oeismatchesbynum;
  cout<<"Analyzing first batch for n up to "<<maxpermsize<<"..."<<endl;
  analyzefirstsequenceset(OEIS, sequencesfilename, setsfilename2, minpermsize, minrevisedsetsize, oeismatchesbynum);
  cout<<"Building sequences for filtered sets up to "<<revisedmaxpermsize<<"..."<<endl;
  if (!data_built) setstosequences(setsfilename2, sequencesfilename2, revisedmaxpermsize);
  unordered_map<int, string> oeismatchesbynum2;
  cout<<"Analyzing filtered sequences"<<endl;
  cout<<"Analyzing second batch for n up to "<<maxpermsize<<"..."<<endl;
  analyzesecondsequenceset(OEIS, sequencesfilename2, minpermsize, oeismatchesbynum2);
  cout<<"Rebuilding OEIS for new-sized sequences..."<<endl;
  Oeis OEIS2("/data/williamkuszmaul/stripped", "/data/williamkuszmaul/names", revisedmaxpermsize - minpermsize + 1, 15); // Note: we allow sequences to start in any of positions 1, 2, ..., 15
  cout<<"Analyzing second batch for n up to "<<revisedmaxpermsize<<"..."<<endl;
  unordered_map<int, string> oeismatchesbynum3;
  analyzesecondsequenceset(OEIS2, sequencesfilename2, minpermsize, oeismatchesbynum3);
  compareoeismatches(oeismatchesbynum2, oeismatchesbynum3);
  cout<<"Writing to output files..."<<endl;
  displaymatches(OEIS2, oeismatchesbynum3, revisedall, true);
  displaymatches(OEIS2, oeismatchesbynum3, revisedbrief, false);
  return 0;  
}
