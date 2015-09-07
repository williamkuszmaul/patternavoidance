// A simple program for speed testing countpatterns.h and fastavoidance.h

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
#include "countpatterns.h"
#include "perm.h"
using namespace std;

// given a file, speed-tests countpatterns on that file
void countallfromfiletest(ifstream &infile, int maxpermsize, bool verbose) {
  string line;
  while (getline(infile, line)) {
    if (verbose) cout<<line<<endl;
    vector < uint64_t > numavoiders;
    timestamp_t start_time = get_timestamp();
    vector < vector < int > > tally;
    vector < vector < int > > completelist;
    countpatterns(line, maxpermsize, tally, completelist, verbose, true);
    timestamp_t end_time = get_timestamp();
    cout<<tally[maxpermsize][0]<<" "<<tally[maxpermsize][1]<<endl;
    //cout<<completelist[maxpermsize][0]<<endl;
    if (verbose) cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  }
  return;
}


// arguments can either be:
// (1) av <input file name> <output file name> <permsize> <whether or not to be verbose>
// (2) cnt <input file name> <permsize> <whether to be verbose> 
int main(int argc, char* argv[]) {
  string choice = argv[1];
  if (choice == "av") { // in this case, speed-test countaviodersfromfile
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
    string infile = argv[2];
    int permsize = stoi(argv[3]);
    ifstream input(infile);
    string verbose = argv[4];
    countallfromfiletest(input, permsize, (verbose == "1"));
    cout<<"Complete"<<endl;
  }
  return 0;
}
