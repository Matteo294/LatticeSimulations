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

    vector<double> m2 {0.173913, -0.307692, -0.571429};
    vector<double> g {2.26843, 1.77515, 1.53061};

    Lattice* latt = new Lattice(16, 16); // 16 time points, 16 space points
    latt->latticeInfo();

    for(int i=0; i<3; i++){
        cout << endl << "********** Simulation " << i+1 << " **********" << endl;
        Phi4* sys = new Phi4(latt, m2[i], g[i]);
        sys->systemInfo();
        Simulator* HMC = new Simulator(sys, latt);
        HMC->runMC(NMC);
    }

    return 0;
}