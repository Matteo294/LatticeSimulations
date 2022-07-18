#include "HMCsimulator.h"

HMCsimulator::HMCsimulator(int nsteps, int ntherm, double dt, double T_MD, SU3YangMills* links) : 
rng_ud(dev_ud()), rng_nd(dev_nd()), ud(0., 1.), nd(0.0, 1.0),
pi (links->Nt, std::vector<std::vector<mat>> (links->Nx, std::vector<mat> (2)))
{this->nsteps=nsteps; this->ntherm=ntherm; this->dt=dt; this->T_MD=T_MD; this->links=links; this->accepted=0;}

HMCsimulator::~HMCsimulator() {;}

double HMCsimulator::MCstep(){

    generateMomenta();

    double deltaH = 0.;
    double oldH = computeHamiltonian();
    auto UCopy = links->U;

    // Evolve with molecular dynamics
    // First step
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            pi[nt][nx][0] += 0.5*dt*evaluateMDdrift(nt, nx, 0); // along t
            pi[nt][nx][1] += 0.5*dt*evaluateMDdrift(nt, nx, 1); // along x
        }
    }

    // Intermediate steps
    for(double t=dt; t<=T_MD; t+=dt){
        leapfrogStep();
    }

    // Final step
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            links->U[nt][nx][0] += dt*computeUevol(im*dt*pi[nt][nx][0], 20);
            links->U[nt][nx][1] += dt*computeUevol(im*dt*pi[nt][nx][1], 20);

        }
    }
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx = 0; nx<links->Nx; nx++){
            pi[nt][nx][0] += 0.5*dt*evaluateMDdrift(nt, nx, 0);
            pi[nt][nx][1] += 0.5*dt*evaluateMDdrift(nt, nx, 1);
        }
    }

    deltaH = computeHamiltonian() - oldH;
    std::cout << "DeltaH: " << deltaH << std::endl;

    if (deltaH>0){
        double r = ud(dev_ud);
        if (exp(-deltaH) <= r){links->U=UCopy;}
    } else accepted++;
    return deltaH;
}

void HMCsimulator::leapfrogStep(){
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            links->U[nt][nx][0] += computeUevol(im*dt*pi[nt][nx][0], 20);
            links->U[nt][nx][1] += computeUevol(im*dt*pi[nt][nx][1], 20);
        }
    }
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            pi[nt][nx][0] += dt*evaluateMDdrift(nt, nx, 0);
            pi[nt][nx][1] += dt*evaluateMDdrift(nt, nx, 1);
        }
    }
}

// Note: trace is not precisely zero but almost --> need fix?
void HMCsimulator::generateMomenta(){
    double p;
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            pi[nt][nx][0].Zero();
            pi[nt][nx][1].Zero();
            for(int i=0; i<8; i++){
                pi[nt][nx][0] += nd(dev_nd)*links->T[i];
                pi[nt][nx][1] += nd(dev_nd)*links->T[i];
            }
        }
    }
}

double HMCsimulator::computeHamiltonian(){
    std::complex<double> p=0.0, r=0.0;
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            for(int i=0; i<8; i++){
                r = (pi[nt][nx][0]*links->T[i]).trace();
                p += r*r;
                r = (pi[nt][nx][1]*links->T[i]).trace();
                p += r*r;
            }
        }
    }
    return 0.5*p.real() + links->computeAction();
}

mat HMCsimulator::evaluateMDdrift(int nt, int nx, int dir){
    assert(dir==0 || dir==1);
    // Then use directly links->Nt
    auto U = links->U;
    int Nt=links->Nt, Nx=links->Nx;

    mat M;

    if (dir==0) M = U[nt][nx][0] * U[PBCidx(nt+1, Nt)][nx][1] * U[nt][PBCidx(nx+1, Nx)][0].adjoint() * U[nt][nx][1].adjoint();
    else M = U[nt][nx][1] * U[nt][PBCidx(nx+1, Nx)][0] * U[PBCidx(nt+1, Nt)][nx][1].adjoint() * U[nt][nx][0].adjoint();
    return -im*links->beta/12. * (M - M.adjoint());
}

mat HMCsimulator::computeUevol(mat M, int N){return M;}