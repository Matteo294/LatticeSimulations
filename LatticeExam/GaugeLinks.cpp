#include "GaugeLinks.h"

GaugeLinks::GaugeLinks(int Nt, int Nx): U(Nt, std::vector<std::vector<mat>> (Nx, std::vector<mat> (2))){
    this->Nt = Nt;
    this->Nx = Nx;
    for(int nt=0; nt<Nt; nt++){
        for(int nx=0; nx<Nx; nx++){
            for(int dir=0; dir <=1; dir++){
                U[nt][nx][dir] = mat{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}; // Initialize links as identity
            }
        }
    }
}

double GaugeLinks::PolyakovLoop(){
    double P = 0.0;
    mat M {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
    for(int nx=0; nx<Nx; nx++){
        for(int nt=Nt-1; nt>=0; nt--){
            M *= U[nt][nx][0];
        }
        P += M.trace().real();
    }
    return (double) P/(3.0*Nt*Nx);
}

double GaugeLinks::computeAction(){
    double beta=1.0; // USE GLOBAL PARAMETEEEEEEEEER
    int Ncolors=1; // USE GLOBAL PARAMETEEEEEEEEER
    mat plaq;
    double S=0.;
    for(int nx=0; nx<Nx; nx++){
        for(int nt=0; nt<Nt; nt++){
            plaq = U[nt][nx][0] * U[PBCidx(nt+1, Nt)][nx][1] * U[nt][PBCidx(nx+1, Nx)][0].adjoint() * U[nt][nx][1].adjoint();
            S = S - ((double) 1.0/(2.0*Ncolors)*(plaq + plaq.adjoint()).trace()).real(); // Check that S is indeed real
        }
    }
    S = beta*(S + Nt*Nx);
    return S;
}

mat GaugeLinks::evaluateMDdrift(int nt, int nx, int dir){
    assert(dir==0 || dir==1);
    mat D;
    if (dir == 0){
        D = U[PBCidx(nt+1, Nt)][nx][1] * U[nt][PBCidx(PBCidx(nx+1, Nx), Nx)][0].adjoint() * U[nt][nx][1].adjoint();
    } else{
        D = U[PBCidx(nt-1, Nx)][nx][0] * U[PBCidx(nt-1, Nt)][PBCidx(nx+1, Nx)][0].adjoint() * U[PBCidx(nt-1, Nt)][nx][1].adjoint();
    }
    return D;
}