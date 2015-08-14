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
//#include "oeislocal.h"
using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 5) {
    cout<<"Incorrect number of arguments. Needs: size of patterns, number of patterns per set, size of permutation, output file name";
  }
  int patternsize = stoi(argv[1]);
  int numpatterns = stoi(argv[2]);
  int maxpermsize = stoi(argv[3]);
  string sequencesfilename = argv[4];
  ofstream sequencesfile;

  string setsfilename = "patternsetstemp";
  ofstream setsfile;
  setsfile.open(setsfilename, std::ofstream::trunc);
  writepatternsetstofile(setsfile, numpatterns, patternsize, false);
  setsfile.close();

  ifstream setsfilein;
  setsfilein.open(setsfilename);
  sequencesfile.open(sequencesfilename, std::ofstream::trunc);
  countavoidersfromfile(setsfilein, sequencesfile, maxpermsize, false);
  setsfilein.close();
  sequencesfile.close();

  return 0;
}
