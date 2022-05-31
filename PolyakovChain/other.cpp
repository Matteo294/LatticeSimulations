#include "other.h"

using namespace std;

std::random_device dev;
std::mt19937 rng(dev());
std::uniform_real_distribution<> dist(0., 1.);

int PBCidx(int n, int N){
    return (n + N)%N;
}

cx_mat::fixed<2,2> genRandomSU2(){
    cx_mat::fixed<2,2> M(arma::fill::eye);
    cx_mat::fixed<2,2> sigma1 = {{0, cx_double(1)}, {cx_double(1), 0}};
    cx_mat::fixed<2,2> sigma2 = {{0, -im}, {im, 0}};
    cx_mat::fixed<2,2> sigma3 = {{cx_double(1), 0}, {0, cx_double(-1)}};
    vector<double> a {dist(dev), dist(dev), dist(dev)};
    double anorm = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    //M = M*cx_double(cos(anorm)) + cx_double(a[0]/anorm*im)*sigma1 + cx_double(a[1]/anorm*im)*sigma2 + cx_double(a[2]/anorm*im)*sigma3;
    return M;
}