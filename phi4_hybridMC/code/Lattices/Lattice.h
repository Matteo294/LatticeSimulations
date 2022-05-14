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
        std::vector<std::vector<double>> copyConfiguration(); // returns a copy of the current configuration
        void writeConfiguration(std::vector<std::vector<double>> phi); // forces the field into a given configuration
        int Nt, Nx; // Number of points in time and space
        std::vector<std::vector<double>> phi; // scalar field value
    private:
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        std::vector<int> Npoints;
};