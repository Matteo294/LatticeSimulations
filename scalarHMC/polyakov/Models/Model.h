#pragma once
#include "../Lattices/Lattice.h"
#include "../Other/other.h"
#include <vector>

using namespace std;

class Model{
    public:
        Model(class Lattice* l);
        ~Model();
        virtual double evaluateAction(){cout << "parent virtual" << endl;return 0.;}
        virtual void modelInfo(){cout << "parent virtual" << endl;}
        virtual double evaluateMDdrift(int nt, int nx){cout << "parent virtual" << endl; return 0.;}
        virtual double computeHamiltonian(){cout << "parent virtual" << endl; return 0.;}
        virtual void writeSavedConfiguration(){cout << "parent virtual" << endl;;}
        virtual void saveConfiguration(){cout << "parent virtual" << endl;;}
        virtual void newConf(int idx, int copyConf=1){cout << "parent virtual" << endl;;}
        vector<vector<double>> pi; // Conjugate momenta
    protected:
        class Lattice* lattice;
        std::random_device rd_gaussian;
        std::mt19937 seed_gaussian;
        std::normal_distribution<> gaussian;
        std::vector<int> Npoints;
        
};