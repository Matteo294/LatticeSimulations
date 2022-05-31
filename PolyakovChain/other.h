#pragma once 
#include <iostream>
#include <vector>
#include <complex.h>
#include <random>
#include <cmath>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

using namespace std;
using namespace arma;

int PBCidx(int n, int N); // Fix indices for periodic boundary conditions
cx_mat::fixed<2,2> genRandomSU2(); // generate a random SU(2) matrix
const complex<double> im (0,1);