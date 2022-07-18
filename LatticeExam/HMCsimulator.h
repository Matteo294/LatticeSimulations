#pragma once 
#include "SU3YangMills.h"
#include <vector>
#include <random>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <iostream>

using mat = Eigen::Matrix3cd;
using vec3 = Eigen::Vector3d;

class HMCsimulator{
    public:
        HMCsimulator(int nsteps, int ntherm, double dt, double T_MD, SU3YangMills* links); // pass number of MC cycles, number of thermalization steps, integration step molecular dynamics
        ~HMCsimulator(); // destructor

        double MCstep(); // performs one MC step
        void runMC(); // performs the whole simulation (thermalization included)
        void generateMomenta(); // generates a random configuration of momenta for HMC
        void leapfrogStep();
        int PBCidx(int n, int N){return (n+N)%N;}
        double computeHamiltonian();
        double factorial(int n);
        mat evaluateMDdrift(int nt, int nx, int dir);
        mat computeUevol(mat M, int N);

        std::vector<std::vector<std::vector<mat>>> pi; // conjugates momenta for HMC
        SU3YangMills *links; // pointer to object of class SU3YangMills
        int accepted;

        std::random_device dev_ud, dev_nd;
        std::mt19937 rng_ud, rng_nd;
        std::uniform_real_distribution<> ud;
        std::normal_distribution<double> nd;
        
    private:
        int nsteps, ntherm; // number of cycles and number of thermalization cycles
        double dt, T_MD; // integration step and final integration time of molecular dynamics
};