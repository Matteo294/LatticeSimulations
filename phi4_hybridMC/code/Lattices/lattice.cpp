#include "Lattice.h"

using namespace std;

Lattice::~Lattice(){};
Lattice::Lattice(int Nt, int Nx) : Npoints(2, 0){
    this->Nt = Nt;
    this->Nx = Nx;
    this->Npoints[0] = Nt;
    this->Npoints[1] = Nx;
    at = 1;
    ax = 1;
}
void Lattice::setSpacing(double at, double ax){this->at=at; this->ax=ax;}
void Lattice::latticeInfo(){
    cout << "Time points: " << Nt << endl;
    cout << "Space points: " << Nx << endl;
    cout << "Time spacing: " << at << endl;
    cout << "Space spacing: " << ax << endl;
}