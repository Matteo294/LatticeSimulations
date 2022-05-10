#include "Model.h"

using namespace std;

Model::Model(class Lattice* l) : phi(l->Nt, vector<double> (l->Nx, 0.0)), seed_gaussian(rd_gaussian()), gaussian(0, 1){
    this->lattice = l;
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            phi[nt][nx] = gaussian(seed_gaussian);
        }
    }
}
Model::~Model(){};
vector<vector<double>> Model::copyConfiguration(){
    return phi;
}

void Model::writeConfiguration(vector<vector<double>> phi){
    for(int nt = 0; nt < lattice->Nt; nt++){
        for(int nx = 0; nx < lattice->Nx; nx++){
            this->phi[nt][nx] = phi[nt][nx];
        }
    }
}