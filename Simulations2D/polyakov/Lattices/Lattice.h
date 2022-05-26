#pragma once
#include <vector>
#include <iostream>
#include <random>

class Lattice{
    public:
        Lattice(int Nt, int Nx);
        ~Lattice();
        void setSpacing(double at, double ax);
        void latticeInfo();
        int Nt, Nx; // Number of points in time and space 
        double at, ax;       
};