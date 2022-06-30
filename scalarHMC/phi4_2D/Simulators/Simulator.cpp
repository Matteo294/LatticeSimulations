#include "Simulator.h"

using namespace std;

Simulator::Simulator(class Model* s, class Lattice* l) : 

    seed_gaussian(rd_gaussian()), 
    gaussian(0, 1), 
    
    seed_uniform(rd_uniform()),
    uniform(0, 1)
    {acceptance=0; this->lattice=l; this->s=s; this->dt=1e-2; this->T=1.;}

Simulator::~Simulator(){};