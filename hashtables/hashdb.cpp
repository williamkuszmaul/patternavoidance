/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  hashdb.h
 *  A simple hash table implementation
 *
 *  This code is released under the MIT license.
 * 
 *  Copyright (c) 2015 William Kuszmaul
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

#include <vector>
#include "hashdb.h"

static unsigned long long next_power_of_two(unsigned long long v) {
  unsigned long long x = 1;
  while (x < v) x*=2;
  return x;
}

hashdb :: hashdb(unsigned long long startsize)
    : averageinsertiontime(0)
    , array(new perm_t [next_power_of_two(startsize)])
    , maxsize(next_power_of_two(startsize))
    , size(0)
{ //The initial maxsize is startsize
  //cout<<maxsize<<endl;
  memset(array, 0, startsize * sizeof(*array));
}

hashdb :: ~hashdb() {
  delete [] array;
}

unsigned long long hashdb :: getsize () const {
  return size;
}

void hashdb :: add (perm_t perm){
  perm_t element = make_key_nonzero(perm);
  assert_key_nonzero(element);
  if(size >= maxsize/2){ //up to 50% full before we resize
    hashdb temp = hashdb(maxsize*2);
    perm_t *oldarray = array;
    for(unsigned long long x=0; x<maxsize; x++){
      if(not_zero_perm(array[x])) temp.add(revert_stored_key(array[x])); //add all of the original elements
    }
    *this = temp;
    temp.array = oldarray; // so that when the destructor runs on temp, it frees the oldarray, not the new array.  It's ugly, but...
  }
  unsigned long long place = hash_perm(element, maxsize);
  while(not_zero_perm(array[place]) && array[place] != element){
    place=(place+1) & (maxsize - 1);
    averageinsertiontime++;
  }
  if (array[place] != element) {
    size++;
    array[place] = element; // we insert element + 1 into the array
  }
}

bool hashdb :: contains (perm_t perm) const {
  perm_t element = make_key_nonzero(perm);
  assert_key_nonzero(element);
  unsigned long long place = hash_perm(element, maxsize);
  while(not_zero_perm(array[place])){
    if(array[place] == element) return true; //look for element+1 in the array
    place = (place + 1) & (maxsize - 1);
  }
  return false;
}


// unsigned long long hashdb :: hash (unsigned long long key) const
// {
//   key = (~key) + (key << 21); // key = (key << 21) - key - 1;
//   key = key ^ (key >> 24);
//   key = (key + (key << 3)) + (key << 8); // key * 265
//   key = key ^ (key >> 14);
//   key = (key + (key << 2)) + (key << 4); // key * 21
//   key = key ^ (key >> 28);
//   key = key + (key << 31);
//   return key&(maxsize-1); //assume maxsize is power of two
//   //This hash function comes from http://www.concentric.net/~ttwang/tech/inthash.htm
// }

unsigned long long hashdb :: getavtime () const { //just to check how well hash function is doing for data set; should round to zero under normal circumstances
  if(size==0) return 0;
  return averageinsertiontime/size;
}

void hashdb::getvals(vector <perm_t> & vals) const{
  vals.resize(0);
  int count = 0;
  for (int i = 0; i < maxsize; i++) {
    if (not_zero_perm(array[i])) {
      count++;
      vals.push_back(revert_stored_key(array[i]));
    }
  }
  assert(count == size);
}
