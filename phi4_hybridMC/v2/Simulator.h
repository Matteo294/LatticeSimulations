#pragma once
#include <iostream>
#include <random>
#include <vector>
#include "System.h"
#include "Lattice.h"
#include "other.h"

// The hybrid monte carlo could then be implemented as a subclass
class Simulator{
    public:
        Simulator(class System* S, class Lattice* l);
        ~Simulator();
        void runMC(int n);
        void thermalize(int n);
        double computeHamiltonian();
        double Enew, Eold;
        int acceptance;
        class System* s;
        class Lattice* lattice;
        // HMC
        double T; // MD simulation time
        double dt;
        std::vector<std::vector<double>> pi;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        std::random_device rd_uniform;
        std::mt19937 seed_uniform;
        std::uniform_real_distribution<double> uniform;
        void createMomentaFields();
        void leapfrogStep();
        double evaluateDrift(int nt, int nx);

};