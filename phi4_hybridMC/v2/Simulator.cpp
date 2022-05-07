#include "Simulator.h"

using namespace std;

Simulator::Simulator(class System* s, class Lattice* l) : 

    seed_gaussian(rd_gaussian()), 
    gaussian(0, 1), 
    pi(l->Nt, vector<double> (l->Nx, 0.0)),
    seed_uniform(rd_uniform()),
    uniform(0, 1)
    {acceptance=0; this->lattice=l; this->s=s; this->dt=1e-2; this->T=1.;}

Simulator::~Simulator(){};

void Simulator::runMC(int n, double thermalization){
    double deltaE=0., Eproposal=0., E=0., Eold, M;
    double avgdE=0., avgexpdeltaE=0.; // observables to compute
    int Nthermalization = thermalization*n;
    
    int progress = 0;
    acceptance = 0;

    for(int i=0; i<n; i++){

        // Just to print the progress of the simulation
        if (progress != ((int)((double)i/n*10))) {progress = (int)((double)i/n*10); cout << "Progress: " << progress*10 << "%" << endl;}

        // Assign random momenta
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] = gaussian(seed_gaussian);
            }
        }

        E = computeHamiltonian();

        // Create a copy of the field in case the new configuration will be rejected
        vector<vector<double>> phi2 = s->copyConfiguration();
        
        // Evolve with molecular dynamics
        // First step
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] += 0.5*dt*s->evaluateMDdrift(nt, nx);
            }
        }
        // Intermediate steps
        for(double t=dt; t<=T; t+=dt){
            leapfrogStep();
        }
        // Final step
        for(int nt=0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                s->phi[nt][nx] += dt*pi[nt][nx];
            }
        }
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] += 0.5*dt*s->evaluateMDdrift(nt, nx);
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
                    Mcurr += s->phi[nt][nx];
                }
            }
            M += (double) Mcurr/(lattice->Nt*lattice->Nx);
        }
    }
    // Print observables
    cout << endl << "Observables values" << endl << "Acceptance: " << ((double)acceptance/(n-Nthermalization)*100) << "% deltaE: " << (double) avgdE/(n-Nthermalization) << " exp(-deltaE): " << (double) avgexpdeltaE/(n-Nthermalization) << " M: " << M/(n-Nthermalization) << endl;
}

void Simulator::leapfrogStep(){
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            s->phi[nt][nx] += dt*pi[nt][nx];
        }
    }
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            pi[nt][nx] += dt*s->evaluateMDdrift(nt, nx);
        }
    }
}

double Simulator::computeHamiltonian(){
    double sum = 0.0;
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            sum += pi[nt][nx]*pi[nt][nx];
        }
    }
    return 0.5*sum + s->evaluateAction(); // H + S
}