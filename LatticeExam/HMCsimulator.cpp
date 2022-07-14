#include "HMCsimulator.h"

HMCsimulator::HMCsimulator(int nsteps, int ntherm, double dt, double T_MD, GaugeLinks* links) : 
rng(dev()), dist(0, 1.),
pi(links->Nt, std::vector<std::vector<mat>> (links->Nx, std::vector<mat> (2)))
{this->nsteps=nsteps; this->ntherm=ntherm; this->dt=dt; this->T_MD=T_MD; this->links=links;}

HMCsimulator::~HMCsimulator() {;}

void HMCsimulator::MCstep(){
    double deltaS = 0.;
    double oldS = links->computeAction();
    auto UCopy = links->U;
    // Evolve with molecular dynamics
    // First step
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            pi[nt][nx][0] += 0.5*dt*links->evaluateMDdrift(nt, nx, 0); // along t
            pi[nt][nx][1] += 0.5*dt*links->evaluateMDdrift(nt, nx, 1); // along x
            std::cout << pi[nt][nx][0] << std::endl;
        }
    }
    // Intermediate steps
    for(double t=dt; t<=T_MD; t+=dt){
        leapfrogStep();
    }
    // Final step
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            links->U[nt][nx][0] += dt*pi[nt][nx][0];
        }
    }
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx = 0; nx<links->Nx; nx++){
            links->U[nt][nx][0] += 0.5*dt*links->evaluateMDdrift(nt, nx, 0);
            links->U[nt][nx][1] += 0.5*dt*links->evaluateMDdrift(nt, nx, 1);
        }
    }

    deltaS = links->computeAction() - oldS;
    if (deltaS>0){
        double r = dist(dev);
        if (exp(-deltaS) <= r){links->U=UCopy;}
    }
}


void HMCsimulator::leapfrogStep(){
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            links->U[nt][nx][0] += dt*pi[nt][nx][0];
            links->U[nt][nx][1] += dt*pi[nt][nx][1];
        }
    }
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            links->U[nt][nx][0] += dt*links->evaluateMDdrift(nt, nx, 0);
            links->U[nt][nx][1] += dt*links->evaluateMDdrift(nt, nx, 1);
        }
    }
}

void HMCsimulator::generateMomenta(){
    for(int nt=0; nt<links->Nt; nt++){
        for(int nx=0; nx<links->Nx; nx++){
            pi[nt][nx][0] = mat {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
            pi[nt][nx][1] = mat {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        }
    }
}
