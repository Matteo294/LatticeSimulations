#include "Lattice.h"

using namespace std;

Lattice::~Lattice(){};
Lattice::Lattice(int Nt, int Nx, int Ny, int Nz){
    this->Nt = Nt;
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    
}
void Lattice::setSpacing(double at, double ax, double ay, double az){this->at=at; this->ax=ax; this->ay=ay; this->az=az;}
void Lattice::latticeInfo(){
    cout << "Time points: " << Nt << endl;
    cout << "Space points: " << Nx << " " << Ny << " " << Nz << endl;
    cout << "Time spacing: " << at << endl;
    cout << "Space spacing: " << ax << " " << ay << " " << az << endl;
}

