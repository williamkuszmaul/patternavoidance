/*
 *  hashmap.h
 *  A simple hash map implementation
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

#ifndef _HASHMAP_H 
#define _HASHMAP_H // To avoid header being included twice in complilation process.

using namespace std;
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <cassert>

// NOTE: KEYS ARE ADDED AS 1 + THEMSELVES. SO 0 IS NOT A VALID KEY
class hashmap{
public:
  hashmap(unsigned long long startsize, int payloadsize);
  ~hashmap();
  unsigned long long hash(unsigned long long key) const;
  void add(unsigned long long element, void *payload);
  void* getpayload(unsigned long long element) const;
  inline unsigned long long getsize () const {
    return size;
  }
 private:
  unsigned long long maxsize; //empty spots will be initated at 0
  unsigned long long size; // number of elts presents
  void *array;
  unsigned int payloadsize; // size of payload
  unsigned int stepsize; // payload + key size (payloadsize + sizeof(unsigned long long))
  inline unsigned long long getkey(unsigned int index) const {
    return *(unsigned long long*)((char*)array + stepsize * index);
  }
  
  inline void  setkey(unsigned int index, unsigned long long key) const {
    *(unsigned long long*)((char*)array + stepsize * index) = key;
  }
  
  inline void*  getval(unsigned int index) const {
    return (char*)array + stepsize * index + sizeof(unsigned long long);
  }
  
  inline void*  setval(unsigned int index, void *value) const {
    memcpy((char*)array + stepsize * index + sizeof(unsigned long long), value, payloadsize);
  }
  
};

#endif
