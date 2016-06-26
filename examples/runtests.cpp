// A simple program for speed testing countpatterns.h and fastavoidance.h

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
#include "permutilities.h"
#include "buildpatternsets.h"
#include "oeislocal.h"
#include "countpatterns.h"
#include "perm.h"
using namespace std;

// given a file, speed-tests countpatterns on that file
void countallfromfiletest(ifstream &infile, int maxpermsize, bool verbose) {
  string line;
  while (getline(infile, line)) {
    cout<<"----------------"<<endl;
    cout<<line<<endl;
    
    vector < uint64_t > numavoiders;
    timestamp_t start_time = get_timestamp();
    vector < vector < int > > tally;
    vector < vector < int > > completelist;
    countpatterns(line, maxpermsize, tally, completelist, verbose, true);
    timestamp_t end_time = get_timestamp();
    for (int i = 1; i <= maxpermsize; i++) {
      cout<<"Number of permutations in S_"<<i<<" with 0, 1, ... patterns"<<endl;
      for (int j = 0; j < tally[i].size(); j++) {
        cout<<tally[i][j]<<" ";
      }
      cout<<endl;
    }
    if (verbose) cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0<<endl;
  }
  return;
}


// arguments can either be:
// (1) av <input file name> <output file name> <permsize> <whether or not to be verbose>
// (2) cnt <input file name> <permsize> <whether to be verbose> 
int main(int argc, char* argv[]) {
  string choice = argv[1];
  if (choice == "av") { // in this case, speed-test countaviodersfromfile
    assert (argc == 6);
    string infile = argv[2];
    ifstream input(infile);
    string outfile = argv[3];
    ofstream output;
    output.open(outfile, std::ofstream::trunc);
    int permsize = stoi(argv[4]);
    string verbose = argv[5];
    countavoidersfromfile(input, output, permsize, (verbose == "1"));
    cout<<"Complete"<<endl;
  } else if (choice == "cnt") { // in this case speed test countpatterns
    assert (argc == 5);
    string infile = argv[2];
    int permsize = stoi(argv[3]);
    ifstream input(infile);
    string verbose = argv[4];
    countallfromfiletest(input, permsize, (verbose == "1"));
    cout<<"Complete"<<endl;
  } else {
    cout<<"Exiting due to invalid command-line format... \n";
  }
  return 0;
}
