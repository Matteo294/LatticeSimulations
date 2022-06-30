#pragma once
#include "Model.h"

using namespace std;

class Phi4_2d : public Model{
    public:
        Phi4_2d(class Lattice* latt, double m2, double g);
        ~Phi4_2d(){;}
        void modelInfo();
        void saveConfiguration();
        void writeSavedConfiguration();
        double evaluateAction();
        double evaluateMDdrift(int nt, int nx); // drift term for Molecular Dynamics evolution
        double computeHamiltonian();
        vector<vector<double>> phi; // Scalar fields
    protected:
        vector<vector<double>> phicopy;
        double m2, g; // coefficients phi^2 and phi^4

};