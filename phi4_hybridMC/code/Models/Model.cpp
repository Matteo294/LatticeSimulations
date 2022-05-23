#include "Model.h"

using namespace std;

Model::Model(class Lattice* l) : 
    seed_gaussian(rd_gaussian()), gaussian(0, 1), 
    phi(l->Nt, vector<vector<vector<double>>> (l->Nx, vector<vector<double>> (l->Ny, vector<double> (l->Nz, 0)))), 
    pi(l->Nt, vector<vector<vector<double>>> (l->Nx, vector<vector<double>> (l->Ny, vector<double> (l->Nz, 0))))
    {this->lattice = l;}
Model::~Model(){};
