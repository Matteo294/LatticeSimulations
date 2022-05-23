#pragma once
#include "Model.h"

using namespace std;

class Phi4_2d : public Model{
    public:
        Phi4_2d(class Lattice* latt, double m2, double g);
        ~Phi4_2d(){;}
        void modelInfo();
        vector<vector<vector<vector<double>>>> copyConfiguration();
        void writeConfiguration(vector<vector<vector<vector<double>>>> phi);
        double evaluateAction();
        double evaluateMDdrift(int nt, int nx, int ny, int nz); // drift term for Molecular Dynamics evolution
    protected:
        double m2, g; // coefficients phi^2 and phi^4

};