#ifndef _BUILDINPUT_H 
#define _BUILDINPUT_H // To avoid header being included twice in complilation process.

void writepatternsetstofile(ofstream &file, int setsize, int patternsize, bool verbose);

void writepatternsetstofile(ofstream &file, int patternsize, bool verbose);

#endif 
