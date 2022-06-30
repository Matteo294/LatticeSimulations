#pragma once
#include "Simulator.h"

using namespace std;

class HMC: public Simulator{
    public:
        HMC(class Model* s, class Lattice* l);
        ~HMC();
        void runMC(int n, int Nskip=0, double thermalization=0.1);
    protected:
        void leapfrogStep();

};