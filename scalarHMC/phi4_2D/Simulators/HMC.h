#pragma once
#include "Simulator.h"

using namespace std;

class HMC: public Simulator{
    public:
        HMC(class Model* s, class Lattice* l);
        ~HMC();
        double runMC(int n, double thermalization=0.1);
    protected:
        void leapfrogStep();

};