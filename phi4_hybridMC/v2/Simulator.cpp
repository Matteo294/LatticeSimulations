#include "Simulator.h"

using namespace std;

Simulator::Simulator(class System* s, class Lattice* l) : 

    seed_gaussian(rd_gaussian()), 
    gaussian(0, 1), 
    pi(l->Nt, vector<double> (l->Nx, 0.0)),
    seed_uniform(rd_uniform()),
    uniform(0, 1)
    {acceptance=0; this->lattice=l; this->s=s; this->dt=1e-3; this->T=1e-1;}

Simulator::~Simulator(){};

void Simulator::thermalize(int n){
    double deltaE;
    Eold = computeHamiltonian();
    Enew = Eold;
    int i=0;
    cout << "Thermalization started..." << endl;
    do{
        Eold = Enew;
         // Random momenta
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] = gaussian(seed_gaussian);
            }
        }
        vector<vector<double>> phi2 = s->copyConfiguration();
        // Evolve with molecular dynamics
        for(double t=0; t<=T; t+=dt){
            leapfrogStep();
        }
        Enew = this->computeHamiltonian();
        deltaE = Enew-Eold;
        cout << "delta " << deltaE << " " << Enew << " " << Eold << " ";
        if ((deltaE > 0) && (exp(-(deltaE)) < uniform(seed_uniform))){
            s->writeConfiguration(phi2);
            Enew = Eold;
            deltaE = 0.0;
            cout << "refused" << endl;
        } else cout << "accepted" << endl;
        i++;
    } while(i<n);
    cout << "Thermalization completed." << endl;
}

void Simulator::runMC(int n){
    double deltaE=0.;
    double avgdE=0., avgexpdeltaE=0., M=0.; // observables to compute
    
    int progress = 0;
    acceptance = 0;
    for(int i=0; i<n; i++){
        if (progress != ((int)((double)i/n*10))) {progress = (int)((double)i/n*10); cout << "Progress: " << progress*10 << "%" << endl;}

        Eold = Enew;

        // Random momentum
        for(int nt = 0; nt < lattice->Nt; nt++){
            for(int nx = 0; nx < lattice->Nx; nx++){
                pi[nt][nx] = gaussian(seed_gaussian);
            }
        }
        vector<vector<double>> phi2 = s->copyConfiguration();
        
        // Evolve with molecular dynamics
        for(double t=0; t<=lattice->Nt; t+=dt){
            leapfrogStep();
        }

        Enew = this->computeHamiltonian();
        deltaE = Enew-Eold;
        cout << "delta " << deltaE << " " << Enew << " " << Eold << " ";
        if ((deltaE > 0) && (exp(-(deltaE)) < uniform(seed_uniform))){
            s->writeConfiguration(phi2);
            Enew = Eold;
            deltaE = 0.0;
            cout << "refused" << endl;
        }
        else{
            acceptance++;
            cout << "accepted" << endl;
        }

        // Compute observables
        avgdE += deltaE;
        avgexpdeltaE += exp(-deltaE);
        /*double sum=0.0;
        for(int nt=0; nt<lattice->Nt; nt++){
            for(int nx=0; nx<lattice->Nx; nx++){
                sum += s->phi[nt][nx];
            }
        }
        M += (double) 1./(lattice->Nt*lattice->Nx) * sum;*/
    }
    cout << (int) acceptance/n*100 << " " << (double) avgdE/n << " " << (double) avgexpdeltaE/n << " " << M << endl;

    //return vector<double> {Enew, deltaE, expdeltaE, M};
}

void Simulator::leapfrogStep(){
    
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            pi[nt][nx] = pi[nt][nx] + 0.5*dt*evaluateDrift(nt, nx);
            s->phi[nt][nx] += dt*pi[nt][nx];
            pi[nt][nx] = pi[nt][nx] + 0.5*dt*evaluateDrift(nt, nx);
        }
    }
}

double Simulator::evaluateDrift(int nt, int nx){
    double Nt = lattice->Nt;
    double Nx = lattice->Nx;
    return 0.5 * (s->phi[PBCidx(nt+1, Nt)][nx] + s->phi[PBCidx(nt-1, Nt)][nx] + s->phi[nt][PBCidx(nx+1, Nx)] 
        + s->phi[nt][PBCidx(nx-1, Nx)] - 4*s->phi[nt][nx]) - (s->m2*s->phi[nt][nx] 
        - (double) s->g/6.*pow(s->phi[nt][nx], 3));
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