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
        virtual double computeHamiltonian(){return 0.;}
        virtual vector<vector<double>> copyConfiguration(){return vector<vector<double>> (1, vector<double> (1, 0.));}
        virtual void writeConfiguration(vector<vector<double>> phi){;}
        vector<vector<double>> phi; // Scalar fields
        vector<vector<double>> pi; // Conjugate momenta
    protected:
        class Lattice* lattice;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        std::vector<int> Npoints;
        
};