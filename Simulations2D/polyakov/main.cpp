#include <iostream>
#include <vector>
#include <fstream>
#include "Other/other.h"
#include "Models/Polyakov.h"
#include "Lattices/Lattice.h"
#include <omp.h>

using namespace std;

const int ntherm = 1000;
const int NMC = 100000; // Number of Monte Carlo cycles
const double b = 0.5;
const int Nlinks = 10;
const int Nt = 16;
const int Nx = 16;

int main(){
    srand(time(NULL));
    Lattice* latt = new Lattice(Nt, Nx);
    Polyakov* P = new Polyakov(latt, Nlinks, b);
    double S = 0.0;
    double Z = 0.0;
    double avgS = 0.0;
    vector<vector<complex<double>>> Rmat (2, vector<complex<double>> (2, complex<double> (0.0, 0.0)));
    vector<vector<complex<double>>> A = matrixMult(vector<vector<complex<double>>> {{complex<double> (1,0), complex<double> (0,0)}, {complex<double> (0,0), complex<double> (1,0)}},
    vector<vector<complex<double>>> {{complex<double> (1,0), complex<double> (0,0)}, {complex<double> (0,0), complex<double> (1,0)}});
    Rmat = randomSU2();
    cout << compute2X2Determinant(Rmat) << endl;
    for(int k1=0; k1<2; k1++){
                for(int k2=0; k2<2; k2++){
                    cout << Rmat[k1][k2] << " ";
                }
                cout << endl;
            }
    for(int n=0; n<5000; n++){
        for(int i=0; i<Nlinks; i++){
            Rmat = randomSU2();
            P->U[i] = matrixMult(Rmat, P->U[i]);
        }
    }
    for(int n=0; n<NMC; n++){
        for(int i=0; i<sizeof(P->U)/sizeof(P->U[0]); i++){
            Rmat = randomSU2();
            P->U[i] = matrixMult(Rmat, P->U[i]);
        }
        S = P->evaluateAction();
        Z += exp(-S);
        avgS += exp(-S)*S;
    }
    cout << "Results: " << avgS/Z << endl;

    return 0;
}