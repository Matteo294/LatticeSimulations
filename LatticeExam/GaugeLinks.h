#pragma once 
#include <vector>
#include "Eigen/Core"

using mat = Eigen::Matrix3cd;


class GaugeLinks{
    public:
        GaugeLinks(int Nt, int Nx);
        ~GaugeLinks(){;}
        double PolyakovLoop();
        double computeAction();
        mat evaluateMDdrift(int nt, int nx, int direction);
        int PBCidx(int n, int N){return (n+N)%N;}
        std::vector<std::vector<std::vector<mat>>> U; // (posx, posy, dir, i, j)
        int Nt, Nx;
};