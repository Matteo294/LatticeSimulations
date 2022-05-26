#include "Polyakov.h"

Polyakov::Polyakov(class Lattice* latt, int Nlinks, double beta) : 
    Model(latt),
    U (Nlinks, vector<vector<complex<double>>> (2, vector<complex<double>> (2, complex<double> (1.0, 0.0))))
    {this->Nlinks=Nlinks; this->beta = beta;}

void Polyakov::modelInfo(){ 
    cout << "Evaluation of a Polyakov chain" << endl;
    cout << "Nlinks: " << Nlinks << " Beta: " << beta << endl;
}

double Polyakov::evaluateAction(){
    vector<vector<complex<double>>> M (2, vector<complex<double>> (2, complex<double> (2, 0.0)));
    M = U.back();
    for(int i=sizeof(U)/sizeof(U[0])-2; i>=0; i--){
        M = matrixMult(U[i], M);
    }
    return -0.5*beta*trace(M).real();
}