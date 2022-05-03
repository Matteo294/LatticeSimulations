#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

// remember to impose PBC in computeHamiltonian()
// Too much staff for the generator. Could some variables be cancelled?
// Since we plan to have high acceptance, I prefer to update the proposal configuration directly in the state array and eventually modify it back

using namespace std;

std::random_device rd_gaussian;
std::mt19937 seed_gaussian(rd_gaussian());
std::normal_distribution<> gaussian(0, 1);
std::random_device rd_uniform;
std::mt19937 seed_uniform(rd_uniform());
std::uniform_real_distribution<double> uniform(0, 1);

vector<int> Npoints {16, 16}; // Nt, Nx
vector<vector<double>> phi(Npoints[0], vector<double> (Npoints[1], 0.0)); // real field initialised to all zeros
vector<vector<double>> pi(Npoints[0], vector<double> (Npoints[1], 0.0)); // fictious field initialised to all zeros
const double m = 1.0; // mass of the field quanta
const double g = 1.0; // strength of the interaction (phi4)
const double eps = 1e-2; // Molecular dynamics step-size
const double t_MD = 1.0; // Molecular dynamics simulation time
const int NMC = 10000; // Number of Monte Carlo cycles
int acceptance = 0;

double computeLocalAction(double nt, double nx); // Compute the (nt, nx)-point's contribution to the action
double computeHamiltonian(); // Compute the fictious hamiltonian
void leapfrogStep(double dt); // Evolves phi(t,x) in the fictious time tau
vector<double> Mcstep(double Hold); // Performs one Monte Carlo step
void initializeStuff(); // Set up field before starting the simulation

int main(){

    ofstream datafile;
    datafile.open("data.csv");
    datafile << "Enew,Eold" << endl;

    initializeStuff();
    double Hold = computeHamiltonian();
    double val = 0.0;
    vector<double> x;
    cout << "initial energy: " << Hold << endl;
    // Thermalization
    for(int i=0; i<1000; i++){
        x = Mcstep(Hold);
        Hold = x[0];
    }
    // Computation
    for(int i=0; i<NMC; i++){
        x = Mcstep(Hold); 
        datafile << x[0] << "," << Hold << endl;
        Hold = x[0]; // update value of energy
        val += x[1]; // cumulative observable
    }
    cout << "Val: " << (double) val/NMC << endl;
    cout << "Acceptance: " << (double) (acceptance-1000)/NMC << endl; // remove thermalization steps

    return 0;
}

// Compute the (nt, nx)-point's contribution to the action
double computeLocalAction(double nt, double nx){
    if(nt==0 || nx==0 || nt==Npoints[0]-1 || nx==Npoints[1]-1) cout << "cannot evaluate these points, they're fixed by boundary conditions" << endl;
    return -0.5 * phi[nt][nx] * (phi[nt+1][nx] + phi[nt-1][nx] + phi[nt][nx+1] + phi[nt][nx-1] - 4*phi[nt][nx])
    + 0.5*m*m*phi[nt][nx]*phi[nt][nx] + (double) g/24. * pow(phi[nt][nx], 4); 

}

// Compute the fictious hamiltonian
double computeHamiltonian(){
    double H=0.0;
    for(int nt=1; nt<Npoints[0]-1; nt++){
        for(int nx=1; nx<Npoints[1]-1; nx++){
            H += 0.5*pi[nt][nx]*pi[nt][nx] + computeLocalAction(nt, nx);
        }
    }
    return H;
}

// Evolves phi(t,x) in the fictious time tau
void leapfrogStep(double dt){
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[0]; nx++){
            pi[nt][nx] -= 0.5*dt * (m*m*phi[nt][nx] + (double) g/6.*pow(phi[nt][nx], 3));
            phi[nt][nx] += dt*pi[nt][nx];
            pi[nt][nx] -= 0.5*dt * (m*m*phi[nt][nx] + (double) g/6.*pow(phi[nt][nx], 3));
        }
    }
}

// Performs one Monte Carlo step
vector<double> Mcstep(double Hold){
    double Hnew, observable_value;
    // Random momentum
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[1]; nx++){
            pi[nt][nx] = gaussian(seed_gaussian);
        }
    }
    // Evolve with molecular dynamics
    for(double t=0; t<=t_MD; t+=eps){
        leapfrogStep(eps);
    }
    Hnew = computeHamiltonian();
    if ((Hnew > Hold) && (exp(-Hnew)/exp(-Hold) > uniform(seed_uniform))){
        leapfrogStep(-eps);
        Hnew = Hold;
    }
    else{
        acceptance += 1;
    }

    observable_value = exp(Hold-Hnew);

    return vector<double> {Hnew, observable_value};

}

// Set up field before starting the simulation
void initializeStuff(){
    for(int nt=0; nt<Npoints[0]; nt++){
        for(int nx=0; nx<Npoints[1]; nx++){
            pi[nx][nt] = gaussian(seed_gaussian);
            phi[nx][nt] = gaussian(seed_gaussian);
            pi[nx][nt] = 1.0;
            phi[nx][nt] = 0.0;
        }
    }
}