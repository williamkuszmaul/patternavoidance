// make command (includes debugging): g++ -O0 -std=c++11  permutationavoidance.cpp  -g -o permutationavoidance
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
#include <unordered_set>
#include <queue>
using namespace std;

// All permutations start at 0
class Permutation {
private:
public:
  vector<int> perm;
  vector<int> inverse;
  int length;

  Permutation (int lengthtemp) {
    perm.resize(lengthtemp);
    inverse.resize(lengthtemp);
    length = lengthtemp;
  }
  
  ~Permutation() {
  }
  
  void buildfrominverse(int *inversetemp) {
    for (int i = 0; i < length; i++) {
      inverse[i] = inversetemp[i];
      perm[inverse[i]] = i;
    }
  }
  
  void buildfromperm(int* permtemp) {
    for (int i = 0; i < length; i++) {
      perm[i] = permtemp[i];
      inverse[perm[i]] = i;
    }
  }

  int getlength() {
    return length;
  }
  
  int get(int i) {
    return perm[i];
  }
  
  int getinverse(int i) {
    return inverse[i];
  }
  
  void display() {
    for (int i = 0; i < length; i++) {
      cout<<perm[i]<<" ";
    }
    cout<<endl;
  }
};


vector < int > buildinversioncode(Permutation &perm) {
  //perm.display();
  uint64_t seenletters = 0; // bit map of which letters we've seen so far
  int length = perm.getlength();
  vector <int> code(length); //int *code = new int[length];
  for (int i = perm.getlength() - 1; i >= 0; i--) {
    if (perm.get(i) == 0) code[i] = 0; // since shifting by 64 bits is ill-defined for some reason
    else code[i] = __builtin_popcountll(seenletters << (64 - perm.get(i)));
    // std::bitset<64> bitview(seenletters);
    // std::bitset<64> bitview2(seenletters << (64 - perm.get(i)));
    // cout<<code[i]<<" "<<(64 - perm.get(i))<<" "<<bitview<<" "<<bitview2<<" "<<((seenletters << 63)<<1)<<endl;
    seenletters = seenletters | (1 << perm.get(i));
    //cout<<code[i]<<" ";
  }
  //  cout<<endl;
  return code;
}

// bijects permutation of n things to 0, ..., n! -  1. 
uint64_t permtoint(Permutation & perm) {
  uint64_t answer = 0;
  uint64_t index = 1;
  uint64_t factorial = 1;
  vector<int> code = buildinversioncode(perm);
  for (int i = perm.getlength() - 2; i >= 0; i--) {
    answer += factorial * code[i];
    index++;
    factorial *= index;
  }
  return answer;
}

// If you want to include the size of the permutation in your mapping, use this instead of permtoint
uint64_t permtohash(Permutation & perm) {
  return (permtoint(perm) << 8) + perm.getlength();
}


// pass perm by reference so that destructor not run on a copy of it (which would free its actual memory)
void addinpos(Permutation &perm, int pos, int val, Permutation &returnperm) {
  int length = perm.getlength() + 1;
  int newperm[length];
  int posshift = 0;
  for (int i = 0; i < length; i++) {
    if (i == pos) {
      newperm[i] = val;
      posshift = 1;
    } else {
      newperm[i] = perm.get(i - posshift);
      if (newperm[i] >= val) newperm[i]++;
    }
  }
  returnperm.buildfromperm(newperm);
}



// // pass perm by reference so that destructor not run on a copy of it (which would free its actual memory)
// Permutation addinpos(Permutation &perm, int pos, int val) {
//   int length = perm.getlength() + 1;
//   int newperm[length];
//   int posshift = 0;
//   for (int i = 0; i < length; i++) {
//     if (i == pos) {
//       newperm[i] = val;
//       posshift = 1;
//     } else {
//       newperm[i] = perm.get(i - posshift);
//       if (newperm[i] >= val) newperm[i]++;
//     }
//   }
//   Permutation answer(length);
//   answer.buildfromperm(newperm);
//   return answer;
// }


Permutation killpos(Permutation &perm, int pos) {
  int length = perm.getlength() - 1;
  int newinverse[length];
  int posshift = 0;
  for (int i = 0; i < length + 1; i++) {
    if (perm.getinverse(i) == pos) {
      posshift = 1;
    } else {
      newinverse[i - posshift] = perm.getinverse(i);
      if (newinverse[i - posshift] > pos) newinverse[i - posshift]--;
    }
  }
  Permutation answer(length);
  answer.buildfrominverse(newinverse);
  return answer;
}


// assumes all permutations smaller than perm are in avoidset iff they are avoiders
// checks if perm is avoider; maxavoidsize = max size of element in patternset; patternset comprises patterns to avoid.
bool isavoider(Permutation & perm, int maxavoidsize, std::unordered_set<long long> &avoidset, std::unordered_set<long long> &patternset) {
  if (patternset.find(permtohash(perm)) != patternset.end()) { // if is in set of bad patterns
    // cout<<"FAILURE FOR ";
    // perm.display();
    // int *code = buildinversioncode(perm);
    // cout<<code[0]<<code[1]<<code[2]<<endl;
    // cout<<"DUE TO SELF"<<endl;
    // cout<<permtohash(perm)<<endl;
    return false;
  }
  if (perm.getlength() > 1) { // don't deal with permutations of size zero
    for (int i = 0; i < perm.getlength() && i <= maxavoidsize; i++) {
      Permutation shortperm = killpos(perm, i);
      if (avoidset.find(permtohash(shortperm))  == patternset.end()) {
	// cout<<"FAILURE FOR ";
	// perm.display();
	// cout<<"DUE TO ";
	// shortperm.display();
	// cout<<endl;
	return false; // found a subword not avoiding the relavent subset
      }
    }
  }
  //perm.display();
  return true;
}


int countavoiders(unordered_set<long long> patternset, int maxavoidsize, int maxsize) {
  int numavoiders = 0;
  
  std::unordered_set<long long> avoidset;
  int temp[1] = {0};
  Permutation startperm(1);
  startperm.buildfromperm(temp);

  std::queue<Permutation> avoiderstoextend;
  avoiderstoextend.push(startperm);
  avoidset.insert(permtohash(startperm));
  while (avoiderstoextend.empty() == false) {
    Permutation perm = avoiderstoextend.front();
    if (perm.getlength() >= maxsize) {
      break;
    }    
    avoiderstoextend.pop();
    for (int i = 0; i < perm.getlength() + 1; i++) {
      Permutation extendedperm(perm.getlength() + 1);
      addinpos(perm, 0, i, extendedperm);
      if (isavoider(extendedperm, maxavoidsize, avoidset, patternset)) {
      	if (extendedperm.getlength() == maxsize) numavoiders++;
      	avoiderstoextend.push(extendedperm);
	avoidset.insert(permtohash(extendedperm));
      }
    }
  }
  return numavoiders;
}

int numinversions(Permutation perm) {
  int answer = 0;
  //perm.display();
  uint64_t seenletters = 0; // bit map of which letters we've seen so far
  int length = perm.getlength();
  for (int i = perm.getlength() - 1; i >= 0; i--) {
    if (perm.get(i) != 0) answer += __builtin_popcountll(seenletters << (64 - perm.get(i)));
    seenletters = seenletters | (1 << perm.get(i));
  }
  return answer;
}

void fillinversions(unordered_set<long long> patternset, int maxavoidsize, int maxsize) {
  int inversiontable[maxsize][maxsize * maxsize];
  for (int i = 0; i < maxsize; i++) {
    for (int j = 0; j < maxsize * maxsize; j++) {
      inversiontable[i][j] = 0;
    }
  }

  int numavoiders =0 ;
  std::unordered_set<long long> avoidset;
  int temp[1] = {0};
  Permutation startperm(1);
  startperm.buildfromperm(temp);

  std::queue<Permutation> avoiderstoextend;
  avoiderstoextend.push(startperm);
  avoidset.insert(permtohash(startperm));
  while (avoiderstoextend.empty() == false) {
    Permutation perm = avoiderstoextend.front();
    if (perm.getlength() >= maxsize) {
      break;
    }    
    avoiderstoextend.pop();
    for (int i = 0; i < perm.getlength() + 1; i++) {
      Permutation extendedperm(perm.getlength() + 1);
      addinpos(perm, 0, i, extendedperm);
      if (isavoider(extendedperm, maxavoidsize, avoidset, patternset)) {
	inversiontable[extendedperm.getlength() - 1][numinversions(extendedperm)]++;
      	if (extendedperm.getlength() == maxsize) numavoiders++;
      	avoiderstoextend.push(extendedperm);
	avoidset.insert(permtohash(extendedperm));
      }
    }
  }

  for (int i = 1; i < maxsize; i++) {
    for (int j = 0; j < maxsize * (maxsize-1) / 2 + 1; j++) {
      cout<<inversiontable[i][j]<<" ";
      assert(inversiontable[i-1][j] <= inversiontable[i][j]);
    }
    cout<<endl;
  }
  cout<<"--"<<endl;

}

bool isid(Permutation perm) {
  for (int i = 1; i < perm.getlength(); i++) {
    if (perm.get(i-1) > perm.get(i)) return false;
  }
  return true;
}

void buildpermutations(Permutation perm, int currentsize, int finalsize) {
  if (currentsize < finalsize) {
    for (int i = 0; i < currentsize + 1; i++) {
      Permutation extendedperm(perm.getlength() + 1);
      addinpos(perm, 0, i, extendedperm);
      buildpermutations(extendedperm, currentsize + 1, finalsize);
    }
  } else {
    if (perm.get(0) == 0 && perm.get(finalsize - 1) == finalsize - 1 && !isid(perm)) {
      perm.display();
      int n = 11;
      int avoidsize = finalsize;
      std::unordered_set<long long> patternset;
      patternset.insert(permtohash(perm));
      fillinversions(patternset, avoidsize, n);
    }
  }
}

void buildpermutationsstart(int finalsize) {
  int temp[1] = {0};
  Permutation startperm(1);
  startperm.buildfromperm(temp);
  buildpermutations(startperm, 1, finalsize);
}


int main() {
  buildpermutationsstart(4);
  return 0;
  // int n = 10;
  // int avoidsize = 5;
  // int avoid1[] = {0, 1, 3, 2, 4, 5};
  // // int avoid2[] = {3, 1, 2, 0};
  // int maxavoidsize = avoidsize;
  // Permutation permavoid1(avoidsize);
  // //Permutation permavoid2(4);
  // permavoid1.buildfromperm(avoid1);
  // //permavoid2.buildfromperm(avoid2);
  // std::unordered_set<long long> patternset;
  // patternset.insert(permtohash(permavoid1));
  // //patternset.insert(permtohash(permavoid2));
  // //cout<<countavoiders(patternset, maxavoidsize, n)<<endl;
  // fillinversions(patternset, maxavoidsize, n);

  //  return 0;
  // Permutation test(5);
  // int perm[] = {3, 2, 4, 1, 0};
  // test.buildfromperm(perm);
  // buildinversioncode(test);
  // cout<<permtoint(test)<<endl;
  // cout<<permtohash(test)<<endl;
  
  return 0;
}
