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
        std::vector<std::vector<double>> copyConfiguration(); // returns a copy of the current configuration
        void writeConfiguration(std::vector<std::vector<double>> phi); // forces the field into a given configuration
        std::vector<std::vector<double>> phi;
        class Lattice* lattice;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        virtual double evaluateAction(){return 0.;}
        virtual void systemInfo(){};
        virtual double evaluateMDdrift(int nt, int nx){return 0.;}
        
};