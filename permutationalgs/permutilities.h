#ifndef _UTILITIES_H 
#define _UTILITIES_H // To avoid header being included twice in complilation process.

#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>      /* printf, scanf, puts, NULL */
#include <cstdlib>
#include <cstring>
#include <ctime> 
#include <iostream>
#include <fstream>
#include <queue>
#include <sys/time.h>
#include <vector>
#include "perm.h"
#include "hashdb.h"
using namespace std;


// For patterns in S_{<10} can use like this:
//  hashdb patternset(1<<3); 
//  string permlist = "3124 4123 3142 4132";
//  makepatterns(permlist, patternset);
// Fills patternset with permutations in permlist, each as a 64 bit integer
void makepatterns(string permlist, hashdb &patternset, int &maxpatternsize);

// Build prefix table containing all complements of normalizations of prefixes of every perm in permset
void addprefixes(const hashdb &permset, hashdb &table);



// In addition to implementing functions useful for permutations, in
// this file we implement a hack which allows for compilers without
// Cilk to compile the program as long as the make file defines the
// macro _NO_CILK. We do this here since permutilities is included in
// every interesting file.
//
// Additionally, throughout the entire project, all cilk-pragma usages
// and all inclusions of cilk headers have to be surrounded by #ifndef
// _NO_CILK ... #endif in order for this hack to work.
#ifdef _NO_CILK

// Create the cilk namespace.
namespace cilk
{
  // Create the cilk::op_add<T> class.
  template <class T>
    class  op_add {
  public:
    // We define the type val_type to be accessed by reducers
    // containing op_add instances.
    using val_type = T;
    // Remember that op_add elements are always initialized to zero!
    T elt = 0;
    // get_value returns the element by reference.
    inline T& get_value () {
      return elt;
    }
  };

  // Create the cilk::op_list_append class.
  template <class T>
    class  op_list_append {
  public:
    // We define the type val_type to be accessed by reducers
    // containing op_list_append instances.
    using val_type = std::list <T>;
    list <T> elt;
    // get_value returns the element by reference.
    inline list <T>& get_value () {
      return elt;
    }
  };

  // Create the cilk::reducer class.
  template <class  A>
    class reducer {
  public:
    A elt;
    // get_value returns the stored value (not by reference).
    inline typename A::val_type get_value () {
      return elt.get_value ();
    }
    // Dereferencing returns the stored value (by reference).
    inline typename A::val_type& operator *() {
      return elt.get_value ();
    }
  };
}
// We overwrite some cilk primatives:
#define cilk_for for
#define cilk_spawn
#define cilk_sync
#endif



#endif 
