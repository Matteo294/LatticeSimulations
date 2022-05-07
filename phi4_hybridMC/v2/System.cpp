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
double System::evaluateAction(){
    double sum = 0.0;
    int Nt = this->lattice->Nt;
    int Nx = this->lattice->Nx;
    for(int nt=0; nt<Nt; nt++){
        for(int nx=0; nx<Nx; nx++){
            sum = sum - 0.5 * phi[nt][nx] * (phi[PBCidx(nt+1, Nt)][nx] + phi[PBCidx(nt-1, Nt)][nx] + phi[nt][PBCidx(nx+1, Nx)] + phi[nt][PBCidx(nx-1, Nx)] - 4*phi[nt][nx])
            + 0.5*m2*phi[nt][nx]*phi[nt][nx] + (double) g/24. * pow(phi[nt][nx], 4);
        }
    }
    return sum;
}
vector<vector<double>> System::copyConfiguration(){
    return this->phi;
}
void System::writeConfiguration(vector<vector<double>> phi){
    for(int nt = 0; nt < lattice->Nt; nt++){
        for(int nx = 0; nx < lattice->Nx; nx++){
            this->phi[nt][nx] = phi[nt][nx];
        }
    }
}