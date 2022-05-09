#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include "other.h"
#include "Simulator.h"
#include "Models.h"
#include "Lattice.h"

using namespace std;

const int NMC = 10000; // Number of Monte Carlo cycles

int main(){

    ofstream results;
    results.open("results.csv");
    results << "m2,M" << endl;

    double M;

    //vector<double> m2 {0.173913, -0.307692, -0.571429};
    double m2;

    //vector<double> g {2.26843, 1.77515, 1.53061};
    double g = 2.;

    Lattice* latt = new Lattice(16, 16); // 16 time points, 16 space points
    latt->latticeInfo();

    for(int i=0; i<50; i++){
        m2 = (double) (2. - (-2.))/50.*i - 2.; 
        cout << endl << "********** Simulation " << i+1 << " **********" << endl;
        Phi4* sys = new Phi4(latt, m2, g);
        sys->systemInfo();
        Simulator* HMC = new Simulator(sys, latt);
        M = HMC->runMC(NMC);
        results << m2 << "," << M << endl;
    }

    return 0;
}