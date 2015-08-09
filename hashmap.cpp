/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  hashmap.cpp
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

#include "hashmap.h"

static unsigned long long next_power_of_two(unsigned long long v) {
  unsigned long long x = 1;
  while (x < v) x*=2;
  return x;
}

hashmap :: hashmap(unsigned long long startsize, int payloadsize)
  : array(new char [(sizeof(unsigned long long) + payloadsize) * next_power_of_two(startsize)])
  , maxsize(next_power_of_two(startsize))
  , size(0)
  , payloadsize(payloadsize)
  , stepsize(payloadsize + sizeof(unsigned long long))
{
  memset(array, 0, maxsize * stepsize);
}

hashmap :: ~hashmap() {
  delete [] (char*)array;
}


void hashmap :: add (unsigned long long element, void* payload) {
  assert(element != -1L);
  if(size >= maxsize/2){ //up to 50% full before we resize
    hashmap temp = hashmap(maxsize*2, payloadsize);
    void *oldarray = array;
    for(unsigned long long x=0; x<maxsize; x++){
      if(getkey(x) != 0) temp.add(getkey(x), getval(x)); //add all of the original elements
    }
    *this = temp;
    temp.array = oldarray; // so that when the destructor runs on temp, it frees the oldarray, not the new array.  It's ugly, but...
  }
  unsigned long long place = hash(element);
  while(getkey(place)!=0){
    place=(place+1)%maxsize;
  }
  setkey(place, element + 1);
  setval(place, payload);
  size++;
}

void * hashmap :: getpayload(unsigned long long element) const {
    assert(element != -1L);
    unsigned long long place=hash(element);
    while(getkey(place) != 0){
      if(getkey(place) == element+1) return getval(place);
      place=(place+1)%maxsize;
    }
    return NULL;
}

unsigned long long hashmap :: hash (unsigned long long key) const
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key&(maxsize-1); //assume maxsize is power of two
  //This hash function comes from http://www.concentric.net/~ttwang/tech/inthash.htm
}
