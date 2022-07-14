#include <iostream>
#include <complex> 
#include <cmath>
#include "Eigen/Core"
#include "GaugeLinks.h"
#include "HMCsimulator.h"

using namespace std;

const std::complex<double> im (0.0, 1.0);

const double beta=1.0, dt=1e-2, T_MD=1e-1;
const int Npoints=16, Ncolors=1;
const int Nsteps=0, Ntherm=2;

int main(){


    GaugeLinks links(Npoints, Npoints);
    HMCsimulator hmc(Nsteps, Ntherm, dt, T_MD, &links);

    // Thermalization
    for(int n=0; n<Ntherm; n++){
        //cout << links.PolyakovLoop() << endl;
        hmc.generateMomenta();
        hmc.MCstep();
    }

    // Simulation
    double P=0.0, P2=0.0, x;
    for(int n=0; n<Nsteps; n++){
        hmc.generateMomenta();
        hmc.MCstep();
        x = links.PolyakovLoop();
        cout << n << " " << x << endl;
        P += x ;
        P2 += x*x;
    }

    cout << "Polyakov loop: " << P/Nsteps << " +- " << sqrt((P*P/Nsteps/Nsteps - P2/Nsteps)/Nsteps) << endl;

    return 0;
}

/*
for(int i=0; i<links.Nt; i++){
    for(int j=0; j<links.Nx; j++){
        cout << hmc.pi[i][j][0] << endl;
        cout << endl;
        cout << hmc.pi[i][j][1] << endl;
        cout << endl << endl;
    }
}
*/
