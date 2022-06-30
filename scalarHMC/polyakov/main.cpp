#include <iostream>
#include <vector>
#include <fstream>
#include "Other/other.h"
#include "Models/Polyakov.h"
#include "Lattices/Lattice.h"
#include "Simulators/MMC.h"
#include <omp.h>
#include <random>

using namespace std;

const int Nskip = 100;
const double b = 0.9;
const int Nlinks = 10;
const int NMC = 1000000; // Number of Monte Carlo cycles
const int Nt = 16;
const int Nx = 16;

int main(){
    
    srand(time(NULL));
    rand();

    Lattice* latt = new Lattice(Nt, Nx);
    Polyakov* P = new Polyakov(latt, Nlinks, b);
    MMC* sim = new MMC(P, latt);

    sim->runMC(NMC, 100, 0.001);
    

    return 0;
}