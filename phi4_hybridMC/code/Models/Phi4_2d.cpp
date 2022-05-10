#include "Phi4_2d.h"

Phi4_2d::Phi4_2d(class Lattice* latt, double m2, double g) : Model(latt){
    this->m2=m2;
    this->g=g;
}

double Phi4_2d::evaluateAction(){
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

void Phi4_2d::modelInfo(){
    cout << endl << "phi4 theory: -phi*ddphi + m2*phi^2 + g*phi^4" << endl;
    cout << "m2: " << m2 << " g: " << g << endl << endl;
}

double Phi4_2d::evaluateMDdrift(int nt, int nx){
    double Nt = lattice->Nt;
    double Nx = lattice->Nx;
    return phi[PBCidx(nt+1, Nt)][nx] + phi[PBCidx(nt-1, Nt)][nx] + phi[nt][PBCidx(nx+1, Nx)] + phi[nt][PBCidx(nx-1, Nx)] - 4*phi[nt][nx] // kinetic term
        - m2*phi[nt][nx] // phi^2 term
        - (double) g/6.*pow(phi[nt][nx], 3); // phi^4 term
}