#include "Simulator.h"

using namespace std;

Simulator::Simulator(class Model* s, class Lattice* l) : 

    seed_gaussian(rd_gaussian()), 
    gaussian(0, 1), 
    seed_uniform(rd_uniform()),
    uniform(0, 1)
    {acceptance=0; this->lattice=l; this->s=s; this->dt=1e-2; this->T=1.;}

Simulator::~Simulator(){};

double Simulator::runMC(int n, double thermalization){
    double deltaE=0., Eproposal=0., E=0., Eold=0., M=0.;
    double avgdE=0., avgexpdeltaE=0.; // observables to compute
    int Nthermalization = thermalization*n;
    
    //int progress = 0;
    acceptance = 0;

    for(int i=0; i<n; i++){

        // Just to print the progress of the simulation
        //if (progress != ((int)((double)i/n*10))) {progress = (int)((double)i/n*10); cout << "Progress: " << progress*10 << "%" << endl;}

        // Assign random momenta
        for(int nt=0; nt<lattice->Nt; nt++){
            for(int nx=0; nx<lattice->Nx; nx++){
                for(int ny=0; ny<lattice->Ny; ny++){
                    for(int nz=0; nz<lattice->Nz; nz++){
                        s->pi[nt][nx][ny][nz] = gaussian(seed_gaussian);
                    }
                }
            }
        }

        E = computeHamiltonian();

        // Create a copy of the field in case the new configuration will be rejected
        vector<vector<vector<vector<double>>>> phi2 = s->copyConfiguration();
        
        // Evolve with molecular dynamics
        // First step
        for(int nt=0; nt<lattice->Nt; nt++){
            for(int nx=0; nx<lattice->Nx; nx++){
                for(int ny=0; ny<lattice->Ny; ny++){
                    for(int nz=0; nz<lattice->Nz; nz++){
                        s->pi[nt][nx][ny][nz] += 0.5*dt*s->evaluateMDdrift(nt, nx, ny, nz);
                    }
                }
            }
        }
        // Intermediate steps
        for(double t=dt; t<=T; t+=dt){
            leapfrogStep();
        }
        // Final step
        for(int nt=0; nt<lattice->Nt; nt++){
            for(int nx=0; nx<lattice->Nx; nx++){
                for(int ny=0; ny<lattice->Ny; ny++){
                    for(int nz=0; nz<lattice->Nz; nz++){
                        s->phi[nt][nx][ny][nz] += dt*s->pi[nt][nx][ny][nz];
                    }
                }
            }
        }
        for(int nt=0; nt<lattice->Nt; nt++){
            for(int nx=0; nx<lattice->Nx; nx++){
                for(int ny=0; ny<lattice->Ny; ny++){
                    for(int nz=0; nz<lattice->Nz; nz++){
                        s->pi[nt][nx][ny][nz] += 0.5*dt*s->evaluateMDdrift(nt, nx, ny, nz);
                    }
                }
            }
        }
        
        Eproposal = this->computeHamiltonian();
        deltaE = Eproposal-E;

        if ((deltaE > 0) && (exp(-(deltaE)) < uniform(seed_uniform))){
            s->writeConfiguration(phi2);
            deltaE = 0.0;
        }
        else {
            if(i >= Nthermalization) acceptance++;
            E = Eproposal;
        }
        if(i >= Nthermalization){
            // Compute observables
            avgdE += deltaE;
            avgexpdeltaE += exp(-deltaE);
            double Mcurr = 0.;
            for(int nt=0; nt<lattice->Nt; nt++){
                for(int nx=0;nx<lattice->Nx; nx++){
                    for(int ny=0; ny<lattice->Ny; ny++){
                        for(int nz=0; nz<lattice->Nz; nz++){
                            Mcurr += s->phi[nt][nx][ny][nz];
                        }
                    }
                }
            }
            M += (double) Mcurr/(lattice->Nt*lattice->Nx*lattice->Ny*lattice->Nz);
        }
    }
    // Print observables
    return (double)M/(n - Nthermalization);
}

void Simulator::leapfrogStep(){
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            for(int ny=0; ny<lattice->Ny; ny++){
                for(int nz=0; nz<lattice->Nz; nz++){
                    s->phi[nt][nx][ny][nz] += dt*s->pi[nt][nx][ny][nz];
                }
            }
        }
    }
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            for(int ny=0; ny<lattice->Ny; ny++){
                for(int nz=0; nz<lattice->Nz; nz++){
                    s->pi[nt][nx][ny][nz] += dt*s->evaluateMDdrift(nt, nx, ny, nz);
                }
            }
        }
    }
}

double Simulator::computeHamiltonian(){
    double sum = 0.0;
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            for(int ny=0; ny<lattice->Ny; ny++){
                for(int nz=0; nz<lattice->Nz; nz++){
                    sum += s->pi[nt][nx][ny][nz]*s->pi[nt][nx][ny][nz];
                }
            }
        }
    }
    return 0.5*sum + s->evaluateAction(); // H + S
}