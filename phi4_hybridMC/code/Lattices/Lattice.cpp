#include "Lattice.h"

using namespace std;

Lattice::~Lattice(){};
Lattice::Lattice(int Nt, int Nx) : Npoints(2, 0), phi(Nt, vector<double> (Nx, 0.0)), seed_gaussian(rd_gaussian()), gaussian(0, 1){
    this->Nt = Nt;
    this->Nx = Nx;
    this->Npoints[0] = Nt;
    this->Npoints[1] = Nx;
    at = 1;
    ax = 1;
    for(int nt=0; nt<Nt; nt++){
        for(int nx=0; nx<Nx; nx++){
            phi[nt][nx] = gaussian(seed_gaussian);
        }
    }
}
void Lattice::setSpacing(double at, double ax){this->at=at; this->ax=ax;}
void Lattice::latticeInfo(){
    cout << "Time points: " << Nt << endl;
    cout << "Space points: " << Nx << endl;
    cout << "Time spacing: " << at << endl;
    cout << "Space spacing: " << ax << endl;
}
vector<vector<double>> Lattice::copyConfiguration(){
    return phi;
}
void Lattice::writeConfiguration(vector<vector<double>> phi){
    for(int nt = 0; nt < Nt; nt++){
        for(int nx = 0; nx < Nx; nx++){
            this->phi[nt][nx] = phi[nt][nx];
        }
    }
}