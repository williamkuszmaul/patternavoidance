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

int main() {
  int maxpermsize = 10;
  int minpermsize = 5;
  string setsfilename = "testallfours.txt";
  string sequencesfilename = "out-testallfours.txt";
  string oeismatchesfilename = "out-out-testallfours.txt";
  
  ofstream setsfile;
  setsfile.open(setsfilename, std::ofstream::trunc);
  writepatternsetstofile(setsfile, 3, true);
  setsfile.close();

  ifstream setsfilein;
  setsfilein.open(setsfilename);
  ofstream sequencesfile;
  sequencesfile.open(sequencesfilename, std::ofstream::trunc);
  countavoidersfromfile(setsfilein, sequencesfile, maxpermsize, true);
  setsfilein.close();
  sequencesfile.close();

  cout<<"Building local version of OEIS..."<<endl;
  Oeis OEIS("stripped", maxpermsize - minpermsize + 1, 15);
  ifstream sequencesfilein;
  sequencesfilein.open(sequencesfilename);
  ofstream oeismatchesfile;
  oeismatchesfile.open(oeismatchesfilename, std::ofstream::trunc);
  analyzesequencefile(sequencesfilein, oeismatchesfile, minpermsize - 1, OEIS, true);
  sequencesfilein.close();
  oeismatchesfile.close();
}
