#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>


// remember to impose PBC in computeHamiltonian()
// Too much staff for the generator. Could some variables be cancelled?

using namespace std;

std::random_device rd_gaussian;
std::mt19937 seed_gaussian(rd_gaussian());
std::normal_distribution<> gaussian(0, 1);
std::random_device rd_uniform;
std::mt19937 seed_uniform(rd_uniform());
std::uniform_real_distribution<double> uniform(0, 1);

vector<int> Npoints {16, 16}; // Nt, Nx
vector<vector<double>> phi(Npoints[0], vector<double> (Npoints[1], 0.0)); // real field initialised to all zeros
vector<vector<double>> pi(Npoints[0], vector<double> (Npoints[1], 0.0)); // fictious  momenta field initialised to all zeros
vector<vector<double>> phi2(Npoints[0], vector<double> (Npoints[1], 0.0)); // copy of realfield 
vector<vector<double>> pi2(Npoints[0], vector<double> (Npoints[1], 0.0)); // fcopy of fictious momenta field

const double m2 = 0.173913; // mass of the field quanta
const double g = 2.26843; // strength of the interaction (phi4)
const double eps = 0.5e-2; // Molecular dynamics step-size
const double t_MD = 1.0; // Molecular dynamics simulation time
const int NMC = 10000; // Number of Monte Carlo cycles
int acceptance = 0;

double computeActionDensity(double nt, double nx); // Compute the (nt, nx)-point's contribution to the action
double computeHamiltonian(); // Compute the fictious hamiltonian
void leapfrogStep(double dt); // Evolve phi(t,x) in the fictious time tau
vector<double> Mcstep(double Hold); // Performs one Monte Carlo step - returns new value of energy and observables
void initializeStuff(); // Set up field before starting the simulation
int PBCidx(int n, int N); // Fix indices for periodic boundary conditions
void copyConfiguration(int reverse); // // Copy pi,phi into pi2,phi2 or vice reversa if reverse different from 0

int main(){

    ofstream datafile;
    datafile.open("data.csv");
    datafile << "Enew,deltaH,exp(deltaH),M" << endl;

    initializeStuff();
    double Hold = computeHamiltonian();
    double val = 0.0;
    vector<double> x(4, 0.0);
    vector<double> xsum(4, 0.0);

    // Thermalization
    for(int i=0; i<1000; i++){
        x = Mcstep(Hold);
        Hold = x[0];
    }

    acceptance = 0;
    int progress = 0;
    // Computation
    for(int i=0; i<NMC; i++){
        if ((progress) != ((int)((double)i/NMC*10))) {cout << "Progress: " << progress*10 << "%" << endl; progress=(int)((double)i/NMC*10);}
        x = Mcstep(Hold); // Hnew, deltaH, exp(deltaH), M 
        datafile << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << endl;
        Hold = x[0]; // update value of energy
        xsum[0] += x[0];
        xsum[1] += x[1];
        xsum[2] += x[2];
        xsum[3] += x[3];
    }
    cout << "Observables values: " << (double) xsum[0]/NMC << "," << (double) xsum[1]/NMC << "," << (double) xsum[2]/NMC << "," << (double) xsum[3]/NMC << endl;
    cout << "Acceptance: " << (double) acceptance/NMC << endl; // remove thermalization steps

    return 0;
}

// Compute the (nt, nx)-point's contribution to the action
double computeActionDensity(double nt, double nx){
    return -0.5 * phi[nt][nx] * (phi[PBCidx(nt+1, Npoints[0])][nx] + phi[PBCidx(nt-1, Npoints[0])][nx] + phi[nt][PBCidx(nx+1, Npoints[1])] + phi[nt][PBCidx(nx-1, Npoints[1])] - 4*phi[nt][nx])
    + 0.5*m2*phi[nt][nx]*phi[nt][nx] + (double) g/24. * pow(phi[nt][nx], 4); 

}

// Compute the fictious hamiltonian
double computeHamiltonian(){
    double H=0.0;
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[1]; nx++){
            H += 0.5*pi[nt][nx]*pi[nt][nx] + computeActionDensity(nt, nx);
        }
    }
    return H;
}

// Evolves phi(t,x) in the fictious time tau
void leapfrogStep(double dt){
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[0]; nx++){
            pi[nt][nx] = pi[nt][nx] + 0.5*0.5*dt * (phi[PBCidx(nt+1, Npoints[0])][nx] + phi[PBCidx(nt-1, Npoints[0])][nx] + phi[nt][PBCidx(nx+1, Npoints[1])] + phi[nt][PBCidx(nx-1, Npoints[1])] - 4*phi[nt][nx]) - 0.5*dt * (m2*phi[nt][nx] + (double) g/6.*pow(phi[nt][nx], 3)) ;
            phi[nt][nx] += dt*pi[nt][nx];
            pi[nt][nx] = pi[nt][nx] + 0.5*0.5*dt * phi[nt][nx] * (phi[PBCidx(nt+1, Npoints[0])][nx] + phi[PBCidx(nt-1, Npoints[0])][nx] + phi[nt][PBCidx(nx+1, Npoints[1])] + phi[nt][PBCidx(nx-1, Npoints[1])] - 4*phi[nt][nx]) -0.5*dt * (m2*phi[nt][nx] + (double) g/6.*pow(phi[nt][nx], 3));
        }
    }
}

// Performs one Monte Carlo step
vector<double> Mcstep(double Hold){
    double Hnew, deltaH=0., expdeltaH=0., M=0.;
    // Random momentum
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[1]; nx++){
            pi[nt][nx] = gaussian(seed_gaussian);
        }
    }
    copyConfiguration(0);
    // Evolve with molecular dynamics
    for(double t=0; t<=t_MD; t+=eps){
        leapfrogStep(eps);
    }
    Hnew = computeHamiltonian();
    deltaH = Hnew-Hold;
    if ((deltaH > 0) && (exp(-(deltaH)) < uniform(seed_uniform))){
        copyConfiguration(1);
        Hnew = Hold;
        deltaH = 0.0;
    }
    else{
        acceptance += 1;
    }

    // Compute observables
    expdeltaH = exp(-deltaH);
    double s=0.0;
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[1]; nx++){
            s += phi[nt][nx];
        }
    }
    M = (double) 1./(Npoints[0]*Npoints[1]) * s;

    return vector<double> {Hnew, deltaH, expdeltaH, M};

}

// Set up field before starting the simulation
void initializeStuff(){
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[1]; nx++){
            pi[nx][nt] = gaussian(seed_gaussian);
            phi[nx][nt] = gaussian(seed_gaussian);
        }
    }
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[1]; nx++){
            pi[nt][nx] -= 0.5*(phi[PBCidx(nt+1, Npoints[0])][nx] + phi[PBCidx(nt-1, Npoints[0])][nx] + phi[nt][PBCidx(nx+1, Npoints[1])] + phi[nt][PBCidx(nx-1, Npoints[1])] - 4*phi[nt][nx]) - 0.5*eps * (m2*phi[nt][nx] + (double) g/6.*pow(phi[nt][nx], 3)) ;
        }
    }
}

// Fix indices for periodic boundary conditions
int PBCidx(int n, int N){
    if(n == N) n = 0;
    else if (n == -1) n = N-1;
    else if (n>N || (n<-1)) cout << "Invalid nt index";
    return n;
}

// Copy pi,phi into pi2,phi2 or vice reversa if reverse different from 0
void copyConfiguration(int reverse){
    if (reverse == 0){
        for(int nt=0; nt<Npoints[0]; nt++){
            for(int nx=0; nx<Npoints[0]; nx++){
                phi2[nt][nx] = phi[nt][nx];
                pi2[nt][nx] = pi[nt][nx];
            }
        }
    } else {
        for(int nt=0; nt<Npoints[0]; nt++){
            for(int nx=0; nx<Npoints[0]; nx++){
                phi[nt][nx] = phi2[nt][nx];
                pi[nt][nx] = pi2[nt][nx];
            }
        }
    }
}