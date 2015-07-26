//  g++ -O0 -std=c++11  hashmaptest.cpp  hashmap.cpp hashmap.h -g -o hashmaptest

#include <assert.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <time.h> 
#include <stdlib.h>
#include <bitset>
#include <vector>
#include <stdint.h>
// #include <unordered_set>
#include <queue>
#include "hashmap.h"
#include <sys/time.h>
using namespace std;

int main () {
  unsigned long long keys[1000];
  unsigned char vals[1000];
  for (unsigned long long i = 0; i < 1000; i++) {
    keys[i] = i*i;
    vals[i] = (unsigned char)((i * i) % 100);
  }
  hashmap testtable(10000, 1);
  for (int i = 0; i < 1000; i++) {
    testtable.add(keys[i], vals + i);
    char *val = (char*)testtable.getpayload(keys[i]);
    assert(*val == vals[i]);
  }
  for (int i = 0; i < 1000; i++) {
    char *val = (char*)testtable.getpayload(keys[i]);
    assert(*val == vals[i]);
  }
  cout<<"Test passed!"<<endl;
  return 0;
}
