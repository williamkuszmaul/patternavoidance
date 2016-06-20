#ifndef _BUILDINPUT_H 
#define _BUILDINPUT_H // To avoid header being included twice in complilation process.

// Writes to file all pattern sets of size setsize containing patterns
// of size patternsize. HOWERVER, only writes one pattern set per
// trivial wilf-equivalence class (so that the same computation is not
// done multiple times due to symmetry)
void writepatternsetstofile(ofstream &file, int setsize, int patternsize, bool verbose);

// Writes to file all pattern sets containing patterns of size
// patternsize. HOWERVER, only writes one pattern set per trivial
// wilf-equivalence class (so that the same computation is not done
// multiple times due to symmetry). This is implemented to be faster
// than repeated calls to the function-interface above.
// Note: only implemented for patternsize < 5. (is not feasible to use for patternset >= 5 anyway)
void writepatternsetstofile(ofstream &file, int patternsize, bool verbose);

#endif 
