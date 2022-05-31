#include <iostream>
#include <vector>
#include <random>
#include "other.h"

using namespace std;
using namespace arma;

double MCstep();
double computeAction();
double stdDev();

double S, oldS, deltaS;
const int NMC = 1000000;
const int Ntherm = 1000;
const int Nskip = 100;
const int Nlinks = 10;
const double b = 0.9;
int accepted, counter;

cx_cube::fixed<2,2,Nlinks> U;
auto Ucopy = U.slice(0);

vector<double> Svec ((int) NMC/Nskip, 0.0);

random_device rnddev;
mt19937 gen(rnddev());
uniform_real_distribution<> ud(0.0, 1.0); 
int main(){  

    double S;
    for(int i=0; i<Ntherm; i++){
        MCstep();
    }

    int accepted=0;
    for(int i=0; i<NMC; i++){
        S = MCstep();
        if ((i%Nskip) == 0) Svec[(int) i/Nskip] = computeAction();
    }

    cout << "Average action: " << (double) accumulate(Svec.begin(), Svec.end(), 0.)/Svec.size() << " +- " << stdDev()/sqrt(Svec.size()) << endl;
    return 0;
}

double computeAction(){
    cx_mat::fixed<2,2> M = U.slice(0);
    for(int i=1; i<Nlinks; ++i){
        M *= U.slice(i);
    }
    return -0.5*b * trace(M).real();
}


double MCstep(){
    deltaS = 0.;
    for(int nl=0; nl<Nlinks; nl++){
        oldS = computeAction();
        Ucopy = U[nl];
        U.slice(nl) *= genRandomSU2();
        deltaS = computeAction() - oldS;
        if (deltaS>0){
            double r = ud(gen);
            if (exp(-deltaS) <= r){U.slice(nl)=Ucopy;}
        }
    }
    return S;
}

double stdDev() {
	auto const mean = accumulate(Svec.begin(), Svec.end(), 0.0) / Svec.size();
	double sum = 0.0;
	for (auto const e : Svec)
		sum += (e - mean) * (e - mean);
	sum /= Svec.size() - 1;
	return std::sqrt(sum);
}
