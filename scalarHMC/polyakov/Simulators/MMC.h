#pragma once
#include "Simulator.h"

using namespace std;

class MMC: public Simulator{
    public:
        MMC(class Model* s, class Lattice* l);
        ~MMC();
        void runMC(int n, int Nskip=1, double thermalization=0.1);

};