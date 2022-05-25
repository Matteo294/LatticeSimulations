#pragma once
#include "../Lattices/Lattice.h"
#include "../Other/other.h"
#include <vector>

using namespace std;

class Model{
    public:
        Model(class Lattice* l);
        ~Model();
        virtual double evaluateAction(){return 0.;}
        virtual void modelInfo(){;}
        virtual double evaluateMDdrift(int nt, int nx){return 0.;}
        virtual vector<vector<double>> copyConfiguration(){return vector<vector<double>> (1, vector<double> (1, 0.0));} // returns a copy of the configuration of the system (set correct return type in the child if it is not this by reshaping the vector)
        virtual void writeConfiguration(vector<vector<double>>){;} // write a configuration on the system (set argument in child function)
        vector<vector<double>> phi;
    protected:
        class Lattice* lattice;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        std::vector<int> Npoints;
        
};