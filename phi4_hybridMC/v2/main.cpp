#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include "other.h"
#include "Simulator.h"
#include "System.h"
#include "Lattice.h"

using namespace std;

const int Nthermalization = 1000;
const int NMC = 10000; // Number of Monte Carlo cycles

/*
-> Action calculation is implemented in "System.cpp"
-> Hamiltonian and drift force calculation is implemented in "Simulator.cpp" since they depend on the algorithm (HMC)
-> Some parameters are fixed directly in the class constructor 'cause I didn't have time do it properly (e.g. g, m2, dt ...)
-> thermalizae and runMC do the same thing at the moment, BUT THEY ARE DIFFERENT FUNCTIONS (aka remember to edit both)
-> For the leapfrog I referred to the notes but seems there is a sign mistake (for the drift)
!! A possible error is how I imposed PBC (see action and drift calculation)
*/

int main(){

    Lattice* latt = new Lattice(16, 16); // 16 time points, 16 space points
    System* s = new System(latt);

    Simulator* mc = new Simulator(s, latt);

    latt->latticeInfo();

    mc->thermalize(Nthermalization);
    mc->runMC(NMC);

    return 0;
}