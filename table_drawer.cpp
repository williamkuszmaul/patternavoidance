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
#include "countpatterns.h"
using namespace std;

string bigsets[] = {"231",
		      "3421 2314 4231 2431 3241 3412 3142 1342 2413 2341 ",
		      "41352 32514 45231 42153 32451 25314 35421 52314 23154 53421 34215 15342 45213 23415 42315 41253 13542 23145 45321 24153 14253 42351 43521 25143 34512 34152 53142 24531 25134 24351 54231 24135 23541 13452 45123 34251 45132 52431 21453 14523 41532 43512 35124 23514 42531 14352 35412 23451 24513 13524 42513 25413 13425 35241 31542 43152 34125 35142 32541 35214 12453 25431 52341 53241 45312 52413 25341 51342 32415 34521 31452 31524 43251 53412 31425 14532 24315 41523 ",
		      "625143 341625 435216 235146 256143 563412 265413 645231 253416 516432 154632 453162 356412 623514 541263 523614 236154 531264 561243 235416 524613 623415 615342 634215 631524 312564 143562 641532 514263 653142 425361 624153 632451 136524 534162 152634 341265 265314 264153 351246 354621 264135 234651 342165 342651 356214 235461 625413 413652 453261 142563 532641 345162 561432 321564 562134 516324 451623 324516 143625 241653 235164 326145 243561 263514 356421 126453 521364 642153 541632 431562 124635 163452 652413 124653 426315 412653 316245 231456 326451 354612 652341 564321 124563 123564 413526 125634 514632 635412 365412 263154 361245 153642 634512 416532 461532 635214 234165 356241 345621 642351 514362 253146 453621 624315 532164 235641 314625 316254 462135 351462 564213 364251 152643 635124 524631 342615 314265 143526 643251 124536 254613 125643 635142 354261 234561 142536 325146 324651 156234 463215 234156 341562 236415 526341 534216 324156 254361 625341 562314 354216 164523 634125 416235 246153 431526 643152 136425 612453 512643 514623 465321 261435 426351 231564 452136 534126 451362 146532 326415 145326 146523 516243 436251 425631 341526 564123 361425 341652 135642 241635 245631 623154 563214 164352 435612 134265 154623 351642 412635 653421 314256 125364 653241 426153 523416 613524 516342 163524 451263 362541 231546 361524 246351 631452 631425 164532 463152 341256 542613 541623 356142 625431 564231 654231 325461 316542 543261 421635 456213 542316 423561 231465 134652 435621 346251 246513 365421 251643 142635 325614 526134 361254 461235 546321 415623 431625 425136 146352 364512 346152 436521 253641 215463 354162 265134 125463 562143 564132 213564 561342 346521 523641 645312 145236 134562 413625 461325 562413 526314 541362 215634 351264 264351 563421 546231 635241 135426 513624 245163 645132 254631 416253 365142 361542 241536 641253 316524 536142 416523 432615 641352 652314 264513 215364 516423 246315 542163 245316 625134 316452 315624 532461 251436 365241 156243 426531 425613 426513 435126 632415 236541 364215 264531 624351 531462 153462 352164 325641 153264 642513 641523 351624 416325 246135 364521 652431 462531 136254 634251 256431 542631 561324 614532 165342 432516 263541 451326 562431 462351 315462 342561 645213 243615 314562 635421 624513 154263 421536 145362 461523 314652 521463 456231 251634 352614 145623 214635 431652 524136 512634 315246 326154 415236 146235 324165 516234 154362 265143 245613 251346 513462 256341 231645 452316 361452 456312 362415 214563 623451 362451 423516 135264 415632 436152 241356 243516 423165 614352 465132 241563 251364 153426 261345 265341 156342 136452 461253 465213 623145 234615 263415 526431 261543 463125 521634 364125 362514 624135 362154 563124 543612 523164 513642 362145 436125 436215 563241 513264 524163 614253 245361 432561 526413 415326 234516 346125 451632 261453 146325 614523 243165 236514 352461 254136 254163 634521 216453 325416 625314 156324 453216 324561 632514 461352 145632 613452 351426 342516 132564 152463 261534 253164 251463 623541 564312 326514 563142 263145 536241 356124 421653 352416 315426 324615 143652 463521 463512 561423 245136 453612 365124 546213 536214 643521 621453 346512 265431 456132 451236 462153 243651 412536 256314 435261 543621 524316 436512 532614 462315 423651 524361 214536 135624 134526 415362 365214 536412 315642 416352 352641 534261 314526 546132 425316 523146 364152 465123 465312 316425 631542 645123 345612 315264 264315 613542 642531 345261 241365 231654 236145 521643 254316 215643 634152 263451 253614 256134 135246 546123 423156 345126 534612 613425 561234 546312 134625 162453 415263 136245 536124 435162 526143 163542 542361 562341 235614 653412 153624 452631 452163 243156 146253 164253 421563 426135 145263 214653 352146 452361 531426 523461 136542 346215 261354 532416 253461 513426 512463 236451 463251 163425 543162 536421 345216 425163 423615 531624 643512 156432 354126 531642 256413 642315 325164 246531 651342 456321 645321 412563 342156 152364 134256 624531 453126 534621 465231 462513 452613 413562 326541 456123 142653 156423 135462 632541 432651 512364 "};


uint64_t Catalan(uint64_t n) {
  assert (n > 0);
  uint64_t Catalans[] = {1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786, 208012, 742900, 2674440, 9694845, 35357670, 129644790, 477638700, 1767263190, 6564120420, 24466267020, 91482563640, 343059613650, 1289904147324};
  assert(n < sizeof(Catalans));
  return Catalans[n];
}

uint64_t choose(uint64_t tmpn, uint64_t tmpk) {
  uint64_t n = tmpn - tmpk + 1;
  uint64_t k = tmpk + 1;
  // cout<<tmpn<<" "<<tmpk<<endl;
  uint64_t array[n][k];
  for (int i = 0; i < n; i++) array[i][0] = 1;
  for (int j = 0; j < k; j++) array[0][j] = 1;
  for (int i = 1; i < n; i++) {
    for (int j = 1; j < k; j++) {
      array[i][j] = array[i-1][j] + array[i][j-1];
    }
  }
  return array[n-1][k-1];
}

uint64_t factorial(int n) {
  double answer = 1;
  for (int i = 1; i <= n; i++) answer *= i;
  return answer;
}

// used for case of looking for k-subsequences containing a 231-pattern
double allestimatedseqcheck(int k, int n) {
  double bottom = 0;
  for (int i = 1; i <= k - 3; i++) {
    bottom += choose(n - k + i, i);
  }
  bottom += (((double)factorial(k -  1) - (double)1) / factorial(k - 1)) * (double) choose(n - 2, k - 2);
  return bottom;
}

// used for single pattern in S_k
double singleestimatedseqcheck(int k, int n) {
  double bottom = 0;
  for (int i = 1; i <= k - 2; i++) {
    bottom += (((double)1) / factorial(i + 1)) * choose(n - k + i, i);
  }
  return bottom;
}

template <typename T>
void buildtable(int startk, int startn, int endk, int endn, vector < vector < T > > array) {
    cout<<"\\begin{tabular}{l |";
  for (int k = startk; k <= endk; k++) {
    cout<<" l ";
  }
  cout<<"}"<<endl;
  for (int k = startk; k <= endk; k++) {
    cout<<" & ";
    cout<<k;
  }
  cout<<" \\\\ \\hline"<<endl;
  for (int n = startn; n <= endn; n++) {
    cout<<n;
    for (int k = startk; k <= endk; k++) {
      cout<<" & "<<array[k][n];
    }
    cout<<" \\\\"<<endl;
  }
  cout<<"\\end{tabular}"<<endl;
}


void test231extensions_avoid(bool getstats) {
    int startk = 3; // just start index. will actually start in first element of bitsets
    int endk = 6; 
    int startn = 8;
    int endn = 13;
    vector < vector < double  > > runtimes (endk + 1, vector < double > (endn + 1, 0));
    vector < vector < uint64_t  > > stat1 (endk + 1, vector < uint64_t > (endn + 1));
    vector < vector < uint64_t  > > stat2 (endk + 1, vector < uint64_t > (endn + 1));
    vector < vector < uint64_t  > > stat3 (endk + 1, vector < uint64_t > (endn + 1));
    vector < vector < uint64_t  > > stat4 (endk + 1, vector < uint64_t > (endn + 1));
    
    vector < vector < double > > workperavoider  (endk + 1, vector < double > (endn + 1));
    vector < vector < double > > fractiononavoiders  (endk + 1, vector < double > (endn + 1));
    vector < vector < double > > specialratio  (endk + 1, vector < double > (endn + 1));
    vector < vector < double > > secspercheck  (endk + 1, vector < double > (endn + 1));
    
    for (int n = startn; n <= endn; n++) {
      for (int k = startk; k <= endk; k++) {
	runtimes[k][n] = run_interior_experiment(bigsets[k - startk], n);
	stat1[k][n] = getstat1();
	stat2[k][n] = getstat2();
	stat3[k][n] = getstat3();
	stat4[k][n] = getstat4();
	workperavoider[k][n] = (double)stat2[k][n] / (double)stat3[k][n]; // average number of sequence checks per avoider in S_n (looking ONLY at S_n here)
	fractiononavoiders[k][n] = (double)stat2[k][n] / (double)stat1[k][n]; // percent of sequence checks in S_n spent on avoiders (looking ONLY at S_n here)
	secspercheck[k][n] = ((double)runtimes[k][n] / (double)stat4[k][n]); // average number seconds spent per sequence check in the entire algorithm 
	specialratio[k][n] = workperavoider[k][n] / (singleestimatedseqcheck(k, n) + 1); // ratio of actual checks per avoider in S_n to expected number. The +1 is to include the single check of a subsequence of length 2 counted by stat1
      }
    }
    cout<<"Computations for 231-extensions:"<<endl;
    cout<<"Runtimes: "<<endl;
    buildtable(startk, startn, endk, endn, runtimes);
    if (getstats) {
      cout<<"Percent of sequence checks in S_n spent on avoiders: "<<endl;
      buildtable(startk, startn, endk, endn, fractiononavoiders);
      cout<<"Average number of sequence checks per avoider in S_n:"<<endl;
      buildtable(startk, startn, endk, endn, workperavoider);
      cout<<"Average number of seconds per sequence check (for all permutations):"<<endl;
      buildtable(startk, startn, endk, endn, secspercheck);
      cout<<"Ratio of number of checks per avoider in S_n to number expected: "<<endl;
      buildtable(startk, startn, endk, endn, specialratio);
    }
}

void testsinglepatterns_avoid(bool getstats) {
    string sets[] = {"231", "2431", "24531", "246531"}; // looking at how they compare to equation eqspec, these may be better to use for experimental data
    int startk = 3; // just start index. will actually start in first element of sets
    int endk = 6;
    int startn = 8;
    int endn = 13;
    vector < vector < double  > > runtimes (endk + 1, vector < double > (endn + 1, 0));
    vector < vector < uint64_t  > > stat1 (endk + 1, vector < uint64_t > (endn + 1));
    vector < vector < uint64_t  > > stat2 (endk + 1, vector < uint64_t > (endn + 1));
    vector < vector < uint64_t  > > stat3 (endk + 1, vector < uint64_t > (endn + 1));
    vector < vector < uint64_t  > > stat4 (endk + 1, vector < uint64_t > (endn + 1));
    
    vector < vector < double > > workperavoider  (endk + 1, vector < double > (endn + 1));
    vector < vector < double > > fractiononavoiders  (endk + 1, vector < double > (endn + 1));
    vector < vector < double > > specialratio  (endk + 1, vector < double > (endn + 1));
    vector < vector < double > > secspercheck  (endk + 1, vector < double > (endn + 1));
    
    for (int n = startn; n <= endn; n++) {
      for (int k = startk; k <= endk; k++) {
	runtimes[k][n] = run_interior_experiment(sets[k - startk], n);
	stat1[k][n] = getstat1();
	stat2[k][n] = getstat2();
	stat3[k][n] = getstat3();
	stat4[k][n] = getstat4();
	workperavoider[k][n] = (double)stat2[k][n] / (double)stat3[k][n]; // average number of sequence checks per avoider in S_n (looking ONLY at S_n here)
	fractiononavoiders[k][n] = (double)stat2[k][n] / (double)stat1[k][n]; // percent of sequence checks in S_n spent on avoiders (looking ONLY at S_n here)
	secspercheck[k][n] = ((double)runtimes[k][n] / (double)stat4[k][n]); // average number seconds spent per sequence check in the entire algorithm 
	specialratio[k][n] = workperavoider[k][n] / (singleestimatedseqcheck(k, n) + 1); // ratio of actual checks per avoider in S_n to expected number. The +1 is to include the single check of a subsequence of length 2 counted by stat1
      }
    }
    cout<<"Computations for single patterns:"<<endl;
    cout<<"Runtimes: "<<endl;
    buildtable(startk, startn, endk, endn, runtimes);
    if (getstats) {
      cout<<"Percent of sequence checks in S_n spent on avoiders: "<<endl;
      buildtable(startk, startn, endk, endn, fractiononavoiders);
      cout<<"Average number of sequence checks per avoider in S_n:"<<endl;
      buildtable(startk, startn, endk, endn, workperavoider);
      cout<<"Average number of seconds per sequence check (for all permutations):"<<endl;
      buildtable(startk, startn, endk, endn, secspercheck);
      cout<<"Ratio of number of checks per avoider in S_n to number expected: "<<endl;
      buildtable(startk, startn, endk, endn, specialratio);
    }
}

void testsinglepatterns_count() {
  string sets[] = {"231", "2431", "24531", "246531"}; // looking at how they compare to equation eqspec, these may be better to use for experimental data
  int startk = 3; // just the start index. will actually start in first element of set
  int endk = 6; 
  int startn = 8;
  int endn = 12;
  vector < vector < double  > > runtimes (endk + 1, vector < double > (endn + 1, 0));
  for (int n = startn; n <= endn; n++) {
    for (int k = startk; k <= endk; k++) {
      runtimes[k][n] = run_interior_experiment2(sets[k - startk], n);
    }
  }
  cout<<"Computations for single patterns:"<<endl;
  cout<<"Runtimes: "<<endl;
  buildtable(startk, startn, endk, endn, runtimes);
}

void testmultiplepatterns_count() {
  int startk = 3; // just the start index. will actually start in first element of bigsets
  int endk = 5; 
  int startn = 8;
  int endn = 12;
  vector < vector < double  > > runtimes (endk + 1, vector < double > (endn + 1, 0));
  for (int n = startn; n <= endn; n++) {
    for (int k = startk; k <= endk; k++) {
      runtimes[k][n] = run_interior_experiment2(bigsets[k - startk], n);
    }
  }
  cout<<"Computations for multiple patterns:"<<endl;
  cout<<"Runtimes: "<<endl;
  buildtable(startk, startn, endk, endn, runtimes);
}

int main() {
  // Only use true argument for variants of brute force; will give you stats.
  testsinglepatterns_avoid(false);
  test231extensions_avoid(false);
  testsinglepatterns_count();
  testmultiplepatterns_count();
  return 0;
}
