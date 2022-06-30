#include <iostream>
#include <vector>
#include <fstream>
#include "Other/other.h"
#include "Simulators/HMC.h"
#include "Models/Phi4_2d.h"
#include "Lattices/Lattice.h"
#include <omp.h>

using namespace std;

const int NMC = 10000; // Number of Monte Carlo cycles

int main(){

    ofstream results;
    results.open("results.csv");
    results << "m2,M" << endl;

    double M;

    vector<double> m2 {0.173913, -0.307692, -0.571429};
    double g = 2.;

    Lattice* latt = new Lattice(16, 16); // 16 time points, 16 space points
    latt->latticeInfo();

    #pragma omp parallel for
    for(int i=0; i<3; i++){
        double M;
        Phi4_2d* model = new Phi4_2d(latt, m2[i], g);
        HMC* sim = new HMC(model, latt);
        M = sim->runMC(NMC);
        cout << i << " " << M << endl;
        delete model, sim;
    }

    delete latt;

    return 0;
}