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
  : array(new char [(sizeof(perm_t) + payloadsize) * next_power_of_two(startsize)])
  , maxsize(next_power_of_two(startsize))
  , size(0)
  , payloadsize(payloadsize)
  , stepsize(payloadsize + sizeof(perm_t))
{
  memset(array, 0, maxsize * stepsize);
}

hashmap :: ~hashmap() {
  delete [] (char*)array;
}


void hashmap :: add (perm_t perm, void* payload) {
  perm_t element = make_key_nonzero(perm);
  assert_key_nonzero(element);
  if(size >= maxsize/2){ //up to 50% full before we resize
    hashmap temp = hashmap(maxsize*2, payloadsize);
    void *oldarray = array;
    for(unsigned long long x=0; x<maxsize; x++){
      if(not_zero_perm(getkey(x))) temp.add(revert_stored_key(getkey(x)), getval(x)); //add all of the original elements
    }
    *this = temp;
    temp.array = oldarray; // so that when the destructor runs on temp, it frees the oldarray, not the new array.  It's ugly, but...
  }
  unsigned long long place = hash_perm(element, maxsize);
  while(not_zero_perm(getkey(place)) && getkey(place) != element){
    place=(place+1) & (maxsize - 1);
  }
  if (getkey(place) != element) {
    setkey(place, element);
    setval(place, payload);
    size++;
  }
}

void * hashmap :: getpayload(perm_t perm) const {
  perm_t element = make_key_nonzero(perm);
  assert_key_nonzero(element);
  unsigned long long place=hash_perm(element, maxsize);
  while(not_zero_perm(getkey(place))){
    if(getkey(place) == element) return getval(place);
    place = (place+1) & (maxsize - 1);
  }
  return NULL;
}
