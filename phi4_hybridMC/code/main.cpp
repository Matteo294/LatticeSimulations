#include <iostream>
#include <vector>
#include <fstream>
#include "Other/other.h"
#include "Simulators/Simulator.h"
#include "Models/Model.h"
#include "Models/Phi4_2d.h"
#include "Lattices/Lattice.h"

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
        Phi4_2d* model = new Phi4_2d(latt, m2, g);
        model->modelInfo();
        Simulator* HMC = new Simulator(model, latt);
        M = HMC->runMC(NMC);
        results << m2 << "," << M << endl;
    }

    return 0;
}