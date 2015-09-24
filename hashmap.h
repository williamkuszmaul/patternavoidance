/*
 *  hashmap.h
 *  A very simple hash map implementation for storing pairs
 *  (unsigned long long, payloadsize byte object) Only has insert and
 *  lookup ability. No update or delete.  Uses linear probing,
 *  which often yields better performance than c++ unordered map.
 *  Will throw error if -1 is inserted as a key. (Is only screw case)
 *  If a key is inserted multiple times, it will be stored multiple
 *  times; Is undefined which payload will then be returned for
 *  lookups. (This could easily be made better if it mattered for a
 *  usage case)
 *  Note: Does not call any kind of destructor on payload object when hash table is destructed. So payload object cannot have pointers which need to be freed.
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
#include "perm.h"

// NOTE: KEYS ARE ADDED AS 1 + THEMSELVES. SO -1 IS NOT A VALID KEY
class hashmap{
public:
  hashmap(unsigned long long startsize, int payloadsize);
  ~hashmap();
  void add(perm_t perm, void *payload);
  void* getpayload(perm_t perm) const;
  inline unsigned long long getsize () const {
    return size;
  }
 private:
  unsigned long long maxsize; //empty spots will be initated at 0
  unsigned long long size; // number of elts presents
  void *array;
  unsigned int payloadsize; // size of payload
  unsigned int stepsize; // payload + key size (payloadsize + sizeof(perm_t))
  inline perm_t getkey(unsigned int index) const {
    return *(perm_t*)((char*)array + stepsize * index);
  }
  
  inline void  setkey(unsigned int index, perm_t key) const {
    *(perm_t*)((char*)array + stepsize * index) = key;
  }
  
  inline void*  getval(unsigned int index) const {
    return (char*)array + stepsize * index + sizeof(perm_t);
  }
  
  inline void  setval(unsigned int index, void *value) const {
    memcpy((char*)array + stepsize * index + sizeof(perm_t), value, payloadsize);
  }
  
};

#endif
