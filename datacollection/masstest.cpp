#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
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
#include "countavoiders.h"
#include "buildpatternsets.h"
#include "perm.h"
#include "oeislocal.h"
using namespace std;

void setstosequences(string inputfile, string outputfile, int maxpermsize) {
  ifstream setsfilein;
  setsfilein.open(inputfile);
  ofstream sequencesfile;
  sequencesfile.open(outputfile, std::ofstream::trunc);
  countavoidersfromfile_parallel(setsfilein, sequencesfile, maxpermsize, false);
  setsfilein.close();
  sequencesfile.close();
}

void buildwilfeqsets(string setsfilename, int patternsize, int maxsetsize, int minsetsize) {
  cout<<"Building sets of "<<patternsize<<"-patterns up to trivial Wilf equivalence..."<<endl;
  ofstream setsfile;
  setsfile.open(setsfilename, std::ofstream::trunc);
  for (int setsize = minsetsize; setsize <= maxsetsize; setsize++) {
    writepatternsetstofile(setsfile, setsize, patternsize, false);
  }
  setsfile.close();
}

void analyzefirstsequenceset(Oeis &OEIS, string sequencesfilename, int minpermsize, int minrevisedsetsize, unordered_map<int, string> &oeismatchesbynum) {
  ifstream sequencesfilein;
  sequencesfilein.open(sequencesfilename);
  vector <patternsetinfo> matches;
  int numattempts = 0, numdistinctattempts = 0;
  fillpatternsetinfo(sequencesfilein, OEIS, minpermsize - 1, matches, numattempts, numdistinctattempts); // start with the minpermsize-th entry
  sequencesfilein.close();
  vector <int> numhitswithdeg(5);
  unordered_set<int> hitswithdeg[5];
  for (int i = 0; i < matches.size(); i++) {
    patternsetinfo set = matches[i];
    bool cont = true;
    int degree = 4;
    if (iszeroby(ithderivative(set.sequence, 4), 5)) degree = 3;
    if (iszeroby(ithderivative(set.sequence, 3), 5)) degree = 2;
    if (iszeroby(ithderivative(set.sequence, 2), 5)) degree = 1;
    if (iszeroby(ithderivative(set.sequence, 1), 5)) degree = 0;
    numhitswithdeg[degree]++;
    hitswithdeg[degree].insert(set.oeisnum);
    if (degree == 4) {
      if (oeismatchesbynum.find(set.oeisnum) == oeismatchesbynum.end()) {
	oeismatchesbynum[set.oeisnum] = set.patternset;
      } else {
	oeismatchesbynum[set.oeisnum] = oeismatchesbynum[set.oeisnum] + "\n" + set.patternset;
      }
    }
  }
  cout<<"Total of "<<matches.size()<<" matches out of "<<numattempts<<" attempts, with what appears to be "<<numdistinctattempts<<" total Wilf-classes"<<endl;
  cout<<"A total of "<<hitswithdeg[0].size() + hitswithdeg[1].size() + hitswithdeg[2].size() + hitswithdeg[3].size() + hitswithdeg[4].size()<<" distinct OEIS sequences appear"<<endl;
  cout<<numhitswithdeg[0]<<" matches constant by first derivative position 5"<<endl;
  cout<<numhitswithdeg[1]<<" of remaining matches linear by second derivative position 5"<<endl;
  cout<<numhitswithdeg[2]<<" of remaining matches quadratic by second derivative position 5"<<endl;
  cout<<numhitswithdeg[3]<<" of remaining matches cubic by second derivative position 5"<<endl;
  cout<<"Number distinct distinct oeis sequences matching to constant sequences: "<<hitswithdeg[0].size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over linear sequences: "<<hitswithdeg[1].size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over quadratic sequences: "<<hitswithdeg[2].size()<<endl;
  cout<<"Number distinct distinct oeis sequences matching to left-over cubic sequences: "<<hitswithdeg[3].size()<<endl;
  cout<<"Number distinct oeis sequences for remaining: "<<oeismatchesbynum.size()<<" accounting for "<<numhitswithdeg[4]<<" sequences."<<endl;
  // Note: oeismatchesbynum.size() = hitswithdeg[4].size()
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
     //output<<OEIS.oeisnames[0]<<endl;
     output<<OEIS.oeisnames[iter->first]<<endl;
     output<<numhits<<" times"<<endl;
     if (verbose) output<<iter->second<<endl;
   }

}

uint64_t factorial(int n) {
  double answer = 1;
  for (int i = 1; i <= n; i++) answer *= i;
  return answer;
}


int main() {
  bool data_built = false; // To just run OEIS analysis since sequence files prebuilt
  int minpermsize = 5; // for analysis, not for data collection
  int patternsize = 4;
  int maxpermsize = 13;
  int comparepermsize = 13;
  int minsetsize = 1;
  string setsfilename = "testallfours-sets";
  string sequencesfilename = "testallfours-sequences";
  string revisedall = "testallfours-matches";
  string revisedbrief = "testallfours-matches-brief";
  cout<<"Building sets up to trivial Wilf-equivalence..."<<endl;
  if (!data_built) buildwilfeqsets(setsfilename, patternsize, factorial(patternsize), minsetsize);
  cout<<"Building sequences up to "<<maxpermsize<<"..."<<endl;
  timestamp_t start_time = get_timestamp();
  if (!data_built) setstosequences(setsfilename, sequencesfilename, maxpermsize);
  timestamp_t end_time = get_timestamp();
  if (!data_built) cout<<(end_time - start_time) / 1000000.0<<" seconds for initial computation."<<endl;
  cout<<"Building local version of OEIS..."<<endl;
  // ON OTHER COMPUTERS, STRIPPED AND NAMES FILES WILL HAVE TO BE CORRECTLY REFERRED TO
  Oeis OEIS("stripped", "names", maxpermsize - minpermsize + 1, 15); //Note: we allow sequences to start in any of positions 1, 2, ..., 15
  cout<<"Continuing analysis..."<<endl;
  unordered_map<int, string> oeismatchesbynum;
  cout<<"Analyzing first batch for n up to "<<maxpermsize<<"..."<<endl;
  analyzefirstsequenceset(OEIS, sequencesfilename, minpermsize, minsetsize, oeismatchesbynum);
  cout<<"Building extended local version of OEIS..."<<endl;
  Oeis OEIS2("stripped", "names", comparepermsize - minpermsize + 1, 15); //Note: we allow sequences to start in any of positions 1, 2, ..., 15
  cout<<"Continuing analysis..."<<endl;
  unordered_map<int, string> oeismatchesbynum2;
  cout<<"Analyzing first batch for n up to "<<comparepermsize<<"..."<<endl;
  analyzefirstsequenceset(OEIS2, sequencesfilename, minpermsize, minsetsize, oeismatchesbynum2);
  compareoeismatches(oeismatchesbynum, oeismatchesbynum2);
  cout<<"Writing to output files..."<<endl;
  displaymatches(OEIS2, oeismatchesbynum, revisedall, true);
  displaymatches(OEIS2, oeismatchesbynum, revisedbrief, false);
  return 0;  
}
