#pragma once
#include <vector>
#include <iostream>
#include <random>

class Lattice{
    public:
        Lattice(int Nt, int Nx, int Ny, int Nz);
        ~Lattice();
        void setSpacing(double at, double ax, double ay, double az);
        void latticeInfo();
        int Nt, Nx, Ny, Nz; // Number of points in time and space 
        double at, ax, ay, az;       
};