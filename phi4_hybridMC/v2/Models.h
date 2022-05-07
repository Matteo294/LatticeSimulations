#pragma once
#include "System.h"

using namespace std;

class Phi4 : public System{
    public:
        Phi4(class Lattice* latt, double m2, double g);
        ~Phi4();
        double evaluateAction();
        void systemInfo();
        double m2, g; // coefficients phi^2 and phi^4
        double evaluateMDdrift(int nt, int nx); // drift term for Molecular Dynamics evolution

};