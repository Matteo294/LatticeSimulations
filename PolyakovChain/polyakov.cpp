#include <iostream>
#include <vector>
#include <random>
#include <armadillo>

using namespace std;
using namespace arma;

double MCstep();
double computeAction();
double stdDev();

const   complex<double> im(0.0,1.0);

double S, oldS, deltaS;
const int NMC = 100000;
const int Ntherm = 1000;
const int Nskip = 100;
const int Nlinks = 10;
const double b = 0.9;
int accepted, counter;

cx_cube::fixed<2,2,Nlinks> U;
auto Ucopy = U.slice(0);

vector<double> Svec ((int) NMC/Nskip, 0.0);

std::random_device dev;
std::mt19937 rng(dev());
std::uniform_real_distribution<> dist(0., 1.);

int PBCidx(int n, int N){
    return (n + N)%N;
}

cx_mat::fixed<2,2> genRandomSU2(){
    cx_mat::fixed<2,2> M(arma::fill::eye);
    cx_mat::fixed<2,2> sigma1 = {{0, cx_double(1)}, {cx_double(1), 0}};
    cx_mat::fixed<2,2> sigma2 = {{0, cx_double(0, -1)}, {cx_double(0, 1), 0}};
    cx_mat::fixed<2,2> sigma3 = {{cx_double(1), 0}, {0, cx_double(-1)}};
    vector<double> a {dist(dev), dist(dev), dist(dev)};
    double anorm = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    
    M = cos(anorm)*M + a[0]/anorm*im*sigma1 + a[1]/anorm*im*sigma3;
    //M = cx_double(cos(anorm))*M + cx_double(a[0]/anorm*im)*sigma1 + cx_double(a[1]/anorm*im)*sigma2 + cx_double(a[2]/anorm*im)*sigma3;
    return M;
}

random_device rnddev;
mt19937 gen(rnddev());
uniform_real_distribution<> ud(0.0, 1.0); 
int main(){ 

    for(int i=0; i<Nlinks; i++) U.slice(i) = genRandomSU2();

    double S;
    for(int i=0; i<Ntherm; i++){
        MCstep();
    }

    int accepted=0;
    for(int i=0; i<NMC; i++){
        S = MCstep();
        if ((i%Nskip) == 0) Svec[(int) i/Nskip] = computeAction();
    }

    cout << "Average action: " << (double) accumulate(Svec.begin(), Svec.end(), 0.)/Svec.size() << " +- " << stdDev()/sqrt(Svec.size()) << endl;
    return 0;
}

double computeAction(){
    auto M = U.slice(0);
    for(int i=1; i<Nlinks; ++i){
        M *= U.slice(i);
    }
    return -0.5*b * trace(M).real();
}


double MCstep(){
    deltaS = 0.;
    for(int nl=0; nl<Nlinks; nl++){
        oldS = computeAction();
        Ucopy = U[nl];
        U.slice(nl) *= genRandomSU2();
        deltaS = computeAction() - oldS;
        if (deltaS>0){
            double r = ud(gen);
            if (exp(-deltaS) <= r){U.slice(nl)=Ucopy;}
        }
    }
    return S;
}

double stdDev() {
	auto const mean = accumulate(Svec.begin(), Svec.end(), 0.0) / Svec.size();
	double sum = 0.0;
	for (auto const e : Svec)
		sum += (e - mean) * (e - mean);
	sum /= Svec.size() - 1;
	return sqrt(sum);
}
