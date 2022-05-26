#include "other.h"

using namespace std;

int PBCidx(int n, int N){
    return (n + N)%N;
}

vector<vector<complex<double>>> matrixMult(vector<vector<complex<double>>> A, vector<vector<complex<double>>> B){
    int nrows1 = sizeof(A)/sizeof(A[0]) + 1;
    int ncols1 = sizeof(A[0])/sizeof(A[0][0]) + 1;
    int nrows2 = sizeof(B)/sizeof(B[0]) + 1;
    int ncols2 = sizeof(B[0])/sizeof(B[0][0]) + 1;
    vector<vector<complex<double>>> res(nrows1, vector<complex<double>> (ncols2, 0.0));
    if (ncols1!=nrows2) throw invalid_argument("ncols1 and nrows2 must be equal");
    else{
        for(int i=0; i<nrows1; i++){
            for(int k=0; k<ncols2; k++){
                res[i][k] = 0.0;
                for(int j=0; j<ncols1; j++){
                    res[i][k] += A[i][j]*B[j][k];
                }

            }
        }
    }
    return res;
}

vector<vector<complex<double>>> constMatMult(double a, vector<vector<complex<double>>> A){
    for(int i=0; i<sizeof(A)/sizeof(A[0]); i++){
        for(int j=0; j<sizeof(A[0])/sizeof(A[0][0]); j++){
            A[i][j] *= a;
        }
    }
    return A;
}

vector<vector<complex<double>>> randomSU2(){
    vector<vector<complex<double>>> M(2, vector<complex<double>> (2, 0.0));
    vector<vector<complex<double>>> sigma1 {{0, 1}, {1, 0}};
    vector<vector<complex<double>>> sigma2 {{0, -im}, {im, 0}};
    vector<vector<complex<double>>> sigma3 {{1, 0}, {0, -1}};
    vector<double> a {(double) rand()/RAND_MAX, (double) rand()/RAND_MAX, (double) rand()/RAND_MAX};
    double anorm = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            M[i][j] += im/anorm*sin(anorm)*(a[0]*sigma1[i][j] + a[1]*sigma2[i][j] + a[2]*sigma3[i][j]);
        }
    }
    M[0][0] += cos(anorm);
    M[1][1] += cos(anorm);
    return M;
}

complex<double> trace(vector<vector<complex<double>>> M){
    complex<double> s = complex<double> {0.0, 0.0};
    for(int i=0; i<sizeof(M)/sizeof(M[0]); i++){
            s += M[i][i];
    }
    return s;
}

complex<double> compute2X2Determinant(vector<vector<complex<double>>> M){
    return M[0][0]*M[1][1] - M[1][0]*M[0][1];
}