#pragma once
#include <iostream>
#include <random>
#include <vector>
#include "../Models/Model.h"
#include "../Lattices/Lattice.h"
#include "../Other/other.h"
#include <omp.h>

// The hybrid monte carlo could then be implemented as a subclass
class Simulator{
    public:
        Simulator(class Model* s, class Lattice* l);
        ~Simulator();
        double runMC(int n, double thermalization=0.1);
        double computeHamiltonian();
        int acceptance;
        class Model* s;
        class Lattice* lattice;
        // HMC
        double T; // MD simulation time
        double dt;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        std::random_device rd_uniform;
        std::mt19937 seed_uniform;
        std::uniform_real_distribution<double> uniform;
        void createMomentaFields();
        void leapfrogStep();

};