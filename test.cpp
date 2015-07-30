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
using namespace std;


// for patterns in S_{<10} can use like this:
//  string permlist = "3124 4123 3142 4132";
//  makepatterns(permlist, patternset);
uint64_t makepatterns(string permlist, hashdb &patternset, int &maxpatternsize) {
  maxpatternsize = 0;
  int pos = 0;
  uint64_t perm = 0;
  //cout<<"Pattern set: "<<endl;
  for (int i = 0; i < permlist.size(); i++) {
    if (permlist[i] == ' ') {
      patternset.add(perm);
      //displayperm(perm);
      maxpatternsize = max(pos + 1, maxpatternsize);
      pos = 0;
      perm = 0;
    } else {
      perm = setdigit(perm, pos, (int)(permlist[i] - '0' - 1));
      pos++;
    }
  }
  // No space at end
  patternset.add(perm);
  //displayperm(perm);
  maxpatternsize = max(pos + 1, maxpatternsize);
      
  //cout<<endl;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    cout<<"Program takes argument of input file name and a permutation length"<<endl;
    return 0;
  }

  string filename = argv[1];
  ifstream input;
  input.open(filename);
  vector <string> inputs;
  string line;
  while (getline(input, line)) {
    inputs.push_back(line);
    //cout<<line<<endl;
  }
  input.close();

  string outputfile = "out-"+filename;
  ofstream output;
  output.open(outputfile, std::ofstream::trunc);
  

  int permsize = stoi(argv[2]);
  assert(permsize <= 16);
  for (int i = 0; i < inputs.size(); i++) {
    int maxpatternsize;
    hashdb patternset = hashdb(1<<3);

    //cout<<"Avoid-set: "<<inputs[i]<<endl;
    output<<"#"<<inputs[i]<<endl;
    
    string permlist = inputs[i];
    makepatterns(permlist, patternset, maxpatternsize);
    // Note: could have boolian for whether to fill avoiders vector. And could have boolian for whether to count avoiders of each size. Or just a separte function to count avoiders of each size.
    
    //perm = setdigit(perm, 3, 3);
    timestamp_t start_time = get_timestamp();
    vector < int > numavoiders;
    countavoiders(patternset, maxpatternsize, permsize, numavoiders, (1<<10)); // for large cases, make last argument much larger!
    for (int i = 2; i <= permsize; i++) {
      output<<numavoiders[i]<<" ";
      //cout<<"Number of avoiders of size "<<i<<" is "<<numavoiders[i]<<endl;
    }
    output<<endl;
    //  assert(numavoiders == avoidersvector[permsize].size());
    timestamp_t end_time = get_timestamp();
    //cout<< "Time elapsed (s): "<<(end_time - start_time)/1000000.0L<<endl;
  }
  output.close();
  return 0;
}
