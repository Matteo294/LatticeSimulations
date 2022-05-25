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
        virtual double evaluateMDdrift(int nt, int nx, int ny, int nz){return 0.;}
        virtual vector<vector<vector<vector<double>>>> copyConfiguration(){return vector<vector<vector<vector<double>>>> {{{{0}}}};}
        virtual void writeConfiguration(vector<vector<vector<vector<double>>>>){;} // write a configuration on the system (set argument in child function)
        vector<vector<vector<vector<double>>>> phi;
        vector<vector<vector<vector<double>>>> pi;
    protected:
        class Lattice* lattice;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        std::vector<int> Npoints;
        
};