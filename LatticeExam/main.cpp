#include <iostream>
#include <complex> 
#include <cmath>
#include "Eigen/Core"
#include "SU3YangMills.h"
#include "HMCsimulator.h"
#include "Eigen/Eigenvalues"

using namespace std;

const double b=1.0, dt=1e-2, T_MD=0.1;
const int Npoints=16;
const int Nsteps=500, Ntherm=100;

int main(){


    SU3YangMills model(Npoints, Npoints, b);
    HMCsimulator hmc(Nsteps, Ntherm, dt, T_MD, &model);
    hmc.generateMomenta();

    mat M {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, -2.0}};
    auto eigenv = model.computeEigenvects(M);
    for (const auto& v: eigenv){
        cout << v << endl << endl; 
    }

    /*
    // Thermalization
    for(int n=0; n<Ntherm; n++){
        cout << links.PolyakovLoop() << endl;;
        hmc.MCstep();
    }

    // Simulation
    double P=0.0, P2=0.0, x;
    double avgdeltaH=0.;
    for(int n=0; n<Nsteps; n++){
        avgdeltaH += hmc.MCstep();
        x = links.PolyakovLoop();
        cout << n << " " << x << endl;
        P += x ;
        P2 += x*x;
    }
    
    cout << "Acceptance: " << hmc.accepted << endl;
    cout << "Avg delta H: " << (double) avgdeltaH/Nsteps << endl;
    cout << "Polyakov loop: " << P/Nsteps << " +- " << sqrt((P*P/Nsteps/Nsteps - P2/Nsteps)/Nsteps) << endl;
    */

    return 0;
}

