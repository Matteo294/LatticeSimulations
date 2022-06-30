#pragma once 
#include <iostream>
#include <vector>
#include <complex.h>
#include <random>
#include <cmath>

using namespace std;

int PBCidx(int n, int N); // Fix indices for periodic boundary conditions
vector<vector<complex<double>>> matrixMult(vector<vector<complex<double>>> A, vector<vector<complex<double>>> B); // perform matrix multiplication between A and B
vector<vector<complex<double>>> constMatMult(double a, vector<vector<complex<double>>> A); // multiply matrix by a constant
vector<vector<complex<double>>> randomSU2(); // generate a random SU(2) matrix
const complex<double> im (0,1);
complex<double> trace(vector<vector<complex<double>>> M);
complex<double> compute2X2Determinant(vector<vector<complex<double>>> M);