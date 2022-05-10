#pragma once
#include <vector>
#include <iostream>

class Lattice{
    public:
        Lattice(int Nt, int Nx);
        ~Lattice();
        void setSpacing(double at, double ax);
        void latticeInfo();
        int Nt, Nx; // Number of points in time and space
    private:
        int at, ax; // spacing int time and space
        std::vector<int> Npoints;
};