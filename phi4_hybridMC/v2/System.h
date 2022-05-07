#pragma once
#include <vector>
#include <random>
#include "Lattice.h"
#include "Simulator.h"
#include "other.h"

class System{
    public:
        System(class Lattice* l);
        ~System();
        std::vector<std::vector<double>> copyConfiguration(); // returns a copyt of the current configuration
        void writeConfiguration(std::vector<std::vector<double>> phi); // forces the field into a given configuration
        std::vector<std::vector<double>> phi;
        class Lattice* lattice;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        // This stuff should be put in another class (subclass called phi4)
        double evaluateAction();
        const double m2 = 0.173913; // mass of the field quanta
        const double g = 2.26843; // strength of the interaction (phi4)
        
};