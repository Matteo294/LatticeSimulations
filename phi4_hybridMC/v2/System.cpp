#include "System.h"

using namespace std;

System::System(class Lattice* l) : phi(l->Nt, vector<double> (l->Nx, 0.0)), seed_gaussian(rd_gaussian()), gaussian(0, 1){
    this->lattice = l;
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            phi[nt][nx] = gaussian(seed_gaussian);
        }
    }
}
System::~System(){};
vector<vector<double>> System::copyConfiguration(){
    return phi;
}
void System::writeConfiguration(vector<vector<double>> phi){
    for(int nt = 0; nt < lattice->Nt; nt++){
        for(int nx = 0; nx < lattice->Nx; nx++){
            this->phi[nt][nx] = phi[nt][nx];
        }
    }
}