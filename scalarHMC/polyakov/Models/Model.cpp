#include "Model.h"

using namespace std;

Model::Model(class Lattice* l) : seed_gaussian(rd_gaussian()), gaussian(0, 1), pi(l->Nt, vector<double> (l->Nx, 0.0)) {
    this->lattice = l;
}
Model::~Model(){};
