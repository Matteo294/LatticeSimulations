#include <iostream>
#include <vector>
#include <random>
#include "Eigen/Core"

using mat = Eigen::Matrix2cd;
using vec = Eigen::Vector3d;
using namespace std;

double computeAction();
double MCstep();
double Wilson();

double S, oldS, deltaS;
const int NMC = 1000000, Ntherm = 1000, Nskip = 100;
const int Nt=20, Nx=20;
const double beta=10.;

vector<double> Svec ((int) NMC/Nskip, 0.0);
vector<double> 

random_device rnddev;
mt19937 gen(rnddev());
uniform_real_distribution<> ud(0.0, 1.0);;
mt19937 gen2(rnddev());
uniform_real_distribution<> nd(0.0, 1.0);

int main(){

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

    return 0;
}

double computeAction(){
    return exp(-beta*accumu)
}
double MCstep();
double Wilson();