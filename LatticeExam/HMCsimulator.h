#pragma once 
#include "GaugeLinks.h"
#include <vector>
#include <random>
#include "Eigen/Core"
#include <iostream>

using mat = Eigen::Matrix3cd;

typedef double (*action)(GaugeLinks);

class HMCsimulator{
    public:
        HMCsimulator(int nsteps, int ntherm, double dt, double T_MD, GaugeLinks* links); // pass number of MC cycles, number of thermalization steps, integration step molecular dynamics
        ~HMCsimulator(); // destructor
        void MCstep(); // performs one MC step
        void runMC(); // performs the whole simulation (thermalization included)
        void generateMomenta(); // generates a random configuration of momenta for HMC
        void leapfrogStep();
        int PBCidx(int n, int N){return (n+N)%N;}
        std::vector<std::vector<std::vector<mat>>> pi; // conjugates momenta for HMC
        GaugeLinks *links; // pointer to object of class GaugeLinks
        action S;
        const int seed=10;
        std::random_device dev;
        std::mt19937 rng;
        std::uniform_real_distribution<> dist;
        
    private:
        int nsteps, ntherm; // number of cycles and number of thermalization cycles
        double dt, T_MD; // integration step and final integration time of molecular dynamics
};