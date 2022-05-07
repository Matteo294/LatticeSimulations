#include "other.h"

using namespace std;

// Fix indices for periodic boundary conditions
int PBCidx(int n, int N){
//    if(n == N) n = 0;
//    else if (n == -1) n = N-1;
//    else if (n>N || (n<-1)) cout << "Invalid nt index";
//    return n;

    // alternatively:
    return (n + N)%N;
}
