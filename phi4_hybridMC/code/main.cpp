#include <iostream>
#include <vector>
#include <fstream>
#include "Other/other.h"
#include "Simulators/Simulator.h"
#include "Models/Model.h"
#include "Models/Phi4_2d.h"
#include "Lattices/Lattice.h"
#include <omp.h>

using namespace std;

const int NMC = 16*600; // Number of Monte Carlo cycles

int main(){

    omp_set_num_threads(16);
    int n_per_thread = NMC/16;

    ofstream results;
    results.open("results.csv");
    results << "m2,M" << endl;

    double M;

    vector<double> m2 {0.173913, -0.307692, -0.571429, -0.8};
    double g = 2.;

    Lattice* latt = new Lattice(3, 12, 12, 12); // 16 time points, 16 space points
    latt->latticeInfo();

    for(int i=0; i<4; i++){
        M = 0.;
        #pragma omp parallel reduction(+:M)
        {
            Phi4_2d* model = new Phi4_2d(latt, m2[i], g);
            Simulator* s = new Simulator(model, latt);
            M += s->runMC(n_per_thread);
            delete model, s;
        }
        cout << (double) M/16. << endl;
    }

    delete latt;

    return 0;
}