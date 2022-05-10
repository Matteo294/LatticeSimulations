#pragma once
#include <vector>
#include <random>
#include "../Lattices/Lattice.h"
#include "../Other/other.h"

class Model{
    public:
        Model(class Lattice* l);
        ~Model();
        std::vector<std::vector<double>> copyConfiguration(); // returns a copy of the current configuration
        void writeConfiguration(std::vector<std::vector<double>> phi); // forces the field into a given configuration
        std::vector<std::vector<double>> phi;
        virtual double evaluateAction(){return 0.;}
        virtual void ModelInfo(){};
        virtual double evaluateMDdrift(int nt, int nx){return 0.;}
    protected:
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        class Lattice* lattice;
        
};