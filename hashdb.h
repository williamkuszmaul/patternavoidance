/*
 *  hashdb.h
 *  A simple hash table implementation for storing 64 bit integers.
 *  CANNOT store -1. Will abort if -1 is added or queried.
 *  Uses simple linear probing. In practice, seems to often run much faster than unordered_set<unsigned long long>, which uses chaining.
 *  If same element is added twice, hashdb.cpp will simply store it twice. (although this could easily be changed if needed)
 *
 *  This code is released under the MIT license.
 * 
 *  Copyright (c) 2013 William Kuszmaul
 *  
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *  
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *  
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

#ifndef _HASHDB_H 
#define _HASHDB_H // To avoid header being included twice in complilation process.

using namespace std;
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <cassert>
#include "perm.h"

class hashdb{
public:
  hashdb(unsigned long long startsize);
  ~hashdb();
  void add(perm_t perm);
  bool contains(perm_t perm) const;
  unsigned long long getavtime() const;
  unsigned long long getsize() const;
  // makes vals contain all entries in hash table
  void getvals(vector <perm_t> &vals) const;
 private:
  unsigned long long averageinsertiontime;
  perm_t *array;
  unsigned long long maxsize; //empty spots will be initated at 0 and elements will be added as one plus themselves
  unsigned long long size; 
};

#endif
