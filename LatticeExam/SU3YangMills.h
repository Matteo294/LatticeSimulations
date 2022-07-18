#pragma once 
#include <vector>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <iostream>
#include <random>
#include <cmath>

using mat = Eigen::Matrix3cd;
using vec = Eigen::Vector3cd;

const std::complex<double> im (0.0, 1.0);

class SU3YangMills{
    public:
        SU3YangMills(int Nt, int Nx, double beta);
        ~SU3YangMills(){;}
        double PolyakovLoop();
        double computeAction();
        int PBCidx(int n, int N){return (n+N)%N;}
        std::vector<vec> computeEigenvects(mat M);
        std::vector<std::vector<std::vector<mat>>> U; // (posx, posy, dir, i, j)
        int Nt, Nx;
        double beta;
        const int Ncolors=3;
        std::vector<mat> T; // Gell-Mann matrices
};