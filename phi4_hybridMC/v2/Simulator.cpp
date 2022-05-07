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

void Simulator::thermalize(int n){
    double deltaE;
    Eold = computeHamiltonian();
    Enew = Eold;
    int i=0;
    cout << "Thermalization started..." << endl;
    do{

        // Create a copy of the field in case the new configuration will be rejected
        vector<vector<double>> phi2 = s->copyConfiguration();
    
        // Assign random momenta
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] = gaussian(seed_gaussian);
            }
        }

        Eold = computeHamiltonian();
        
        // Evolve with molecular dynamics
        // First step
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] += 0.5*dt*evaluateDrift(nt, nx);
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
                pi[nt][nx] += 0.5*dt*evaluateDrift(nt, nx);
            }
        }
        
        Enew = this->computeHamiltonian();
        deltaE = Enew-Eold;
        cout << "deltaE " << deltaE << " New energy: " << Enew << " Old energy: " << Eold << " --> ";

        if ((deltaE > 0) && (exp(-(deltaE)) < uniform(seed_uniform))){
            s->writeConfiguration(phi2);
            Enew = Eold;
            deltaE = 0.0;
            cout << "rejected" << endl;
        }
        else{
            acceptance++;
            cout << "accepted" << endl;
        }
        i++;
    } while(i<n);
    cout << "Thermalization completed" << endl;

}

void Simulator::runMC(int n){
    double deltaE=0.;
    double avgdE=0., avgexpdeltaE=0.; // observables to compute
    
    int progress = 0;
    acceptance = 0;

    for(int i=0; i<n; i++){

        // Just to print the progress of the simulation
        if (progress != ((int)((double)i/n*10))) {progress = (int)((double)i/n*10); cout << "Progress: " << progress*10 << "%" << endl;}

        //vector<vector<double>> pi2(pi);
        // Assign random momenta
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] = gaussian(seed_gaussian);
            }
        }

        Eold = computeHamiltonian();

        // Create a copy of the field in case the new configuration will be rejected
        vector<vector<double>> phi2 = s->copyConfiguration();
        
        // Evolve with molecular dynamics
        // First step
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] += 0.5*dt*evaluateDrift(nt, nx);
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
                pi[nt][nx] += 0.5*dt*evaluateDrift(nt, nx);
            }
        }
        
        Enew = this->computeHamiltonian();
        deltaE = Enew-Eold;
        cout << "deltaE " << deltaE << " New energy: " << Enew << " Old energy: " << Eold << " --> ";

        if ((deltaE > 0) && (exp(-(deltaE)) < uniform(seed_uniform))){
            s->writeConfiguration(phi2);
            Enew = Eold;
            deltaE = 0.0;
            cout << "rejected" << endl;
        }
        else{
            acceptance++;
            cout << "accepted" << endl;
        }

        // Compute observables
        avgdE += deltaE;
        avgexpdeltaE += exp(-deltaE);
    }
    // Print observables
    cout << (int) acceptance/n*100 << " " << (double) avgdE/n << " " << (double) avgexpdeltaE/n << endl;
}

void Simulator::leapfrogStep(){
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            s->phi[nt][nx] += dt*pi[nt][nx];
        }
    }
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            pi[nt][nx] += dt*evaluateDrift(nt, nx);
        }
    }
}

double Simulator::evaluateDrift(int nt, int nx){
    double Nt = lattice->Nt;
    double Nx = lattice->Nx;
    return 0;
    return +s->phi[PBCidx(nt+1, Nt)][nx] + s->phi[PBCidx(nt-1, Nt)][nx] + s->phi[nt][PBCidx(nx+1, Nx)] + s->phi[nt][PBCidx(nx-1, Nx)] - 4*s->phi[nt][nx] // kinetic term
        - s->m2*s->phi[nt][nx] // phi^2 term
        - (double) s->g/6.*pow(s->phi[nt][nx], 3); // phi^4 term
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