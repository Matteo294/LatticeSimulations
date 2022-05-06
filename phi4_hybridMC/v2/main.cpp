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
const int NMC = 1000; // Number of Monte Carlo cycles

int main(){

    Lattice* latt = new Lattice(16, 16);
    System* s = new System(latt);

    Simulator* mc = new Simulator(s, latt);

    latt->latticeInfo();

    mc->computeHamiltonian();
    mc->thermalize(Nthermalization);
    mc->runMC(NMC);

    return 0;
}