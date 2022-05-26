#include "Lattice.h"

using namespace std;

Lattice::~Lattice(){};
Lattice::Lattice(int Nt, int Nx){
    this->Nt = Nt;
    this->Nx = Nx;
    
}
void Lattice::setSpacing(double at, double ax){this->at=at; this->ax=ax;}
void Lattice::latticeInfo(){
    cout << "Time points: " << Nt << endl;
    cout << "Space points: " << Nx << endl;
    cout << "Time spacing: " << at << endl;
    cout << "Space spacing: " << ax << endl;
}

