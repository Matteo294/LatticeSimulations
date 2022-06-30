#include "Polyakov.h"

Polyakov::Polyakov(class Lattice* latt, int Nlinks, double beta) : 
    Model(latt),
    U (Nlinks, vector<vector<complex<double>>> (2, vector<complex<double>> (2, complex<double> {1.0, 1.0}))),
    Ucopy (Nlinks, vector<vector<complex<double>>> (2, vector<complex<double>> (2, complex<double> {1.0, 1.0})))
    {this->Nlinks=Nlinks; this->beta = beta;  for(int i=0; i<Nlinks; i++) U[i]=randomSU2();}

void Polyakov::modelInfo(){ 
    cout << "Evaluation of a Polyakov chain" << endl;
    cout << "Nlinks: " << Nlinks << " Beta: " << beta << endl;
}

double Polyakov::evaluateAction(){
    vector<vector<complex<double>>> M (2, vector<complex<double>> (2, complex<double> (2, 0.0)));
    M = U[Nlinks-1];
    for(int i=Nlinks-2; i>=0; i--){
        M = matrixMult(U[i], M);
    }
    return -(double) 0.5*beta*trace(M).real();
}

void Polyakov::printU(int idx){
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            cout << U[idx][i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Polyakov::printUcopy(int idx){
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            cout << Ucopy[idx][i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Polyakov::saveConfiguration(){
    Ucopy = U;
}

void Polyakov::writeSavedConfiguration(){
    U = Ucopy;
}

void Polyakov::newConf(int idx, int copyConf){
    //if (copyConf == 1) saveConfiguration();
    //int l =  (int) (((double) rand()/RAND_MAX) * Nlinks);
    int l = idx; 
    auto Rmat = randomSU2();
    U[l] = matrixMult(Rmat, U[l]);
}