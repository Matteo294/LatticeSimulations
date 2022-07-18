#include "SU3YangMills.h"

SU3YangMills::SU3YangMills(int Nt, int Nx, double beta): 
    U(Nt, std::vector<std::vector<mat>> (Nx, std::vector<mat> (2))), T(8){
    this->Nt=Nt;
    this->Nx=Nx;
    this->beta=beta;
    for(int nt=0; nt<Nt; nt++){
        for(int nx=0; nx<Nx; nx++){
            for(int dir=0; dir <=1; dir++){
                U[nt][nx][dir] = mat{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}; // Initialize links as identity
            }
        }
    }
    /*T[0] = mat {{0., 1., 0.}, {1., 0., 0.}, {0., 0., 0.}};
    T[1] = mat {{0., -im, 0.}, {im, 0., 0.}, {0., 0., 0.}};
    T[2] = mat {{1., 0., 0.}, {0., -1., 0.}, {0., 0., 0.}};
    T[3] = mat {{0., 0., 1.}, {0., 0., 0.}, {1., 0., 0.}};
    T[4] = mat {{0., 0., -im}, {0., 0., 0.}, {im, 0., 0.}};
    T[5] = mat {{0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}};
    T[6] = mat {{0., 0., 0.}, {0., 0., -im}, {0., im, 0.}};
    T[7] = 1./sqrt(3) * mat {{1., 0., 0.}, {0., 1., 0.}, {0., 0., -2.}};*/
}

double SU3YangMills::PolyakovLoop(){
    double P = 0.0;
    mat M {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
    for(int nx=0; nx<Nx; nx++){
        for(int nt=Nt-1; nt>=0; nt--){
            M *= U[nt][nx][0];
        }
        P += M.trace().real();
    }
    return (double) P/(3.0*Nt*Nx);
}

double SU3YangMills::computeAction(){
    mat plaq;
    std::complex<double> r;
    double S=0.;
    for(int nx=0; nx<Nx; nx++){
        for(int nt=0; nt<Nt; nt++){
            plaq = U[nt][nx][0] * U[PBCidx(nt+1, Nt)][nx][1] * U[nt][PBCidx(nx+1, Nx)][0].adjoint() * U[nt][nx][1].adjoint();
            r = 0.;
            for(int i=0; i<8; i++){
                r = ((plaq + plaq.adjoint())).trace();
                S -=  ((double) 1.0/(2.0*Ncolors)*r*r).real(); // Check that S is indeed real
            }
        }
    }
    S = beta*(S + Nt*Nx);
    return S;
}

// add assert to check that the matrix is in su(3)
std::vector<vec> SU3YangMills::computeEigenvects(mat M){
    mat Mcopy = M;
    vec eigenvals;

    // Find eigenvalues
    std::vector<vec> r (3);
    std::vector<vec> v (3);

    // Find eigenvalues solving cubic equation via Cardano's method
    std::complex<double> Q = -(M*M).trace()/6.;
    std::complex<double> R = -M.determinant();
    std::complex<double> S = pow(R + sqrt(Q*Q*Q + R*R), 1./3.);
    std::complex<double> T = pow(R - sqrt(Q*Q*Q + R*R), 1./3.);
    std::complex<double> a = S+T;
    std::complex<double> b = S-T;
    eigenvals = vec{a, -0.5*a + 0.5*im*sqrt(3.)*b, -0.5*a - 0.5*im*sqrt(3.)*b};
    
    // For each eigenvalue calculate the corresponding eigenvector
    for(int eig=0; eig<3; eig++){
        M = Mcopy;
        for(int i=0; i<3; i++){
            M(i,i) = M(i,i) - eigenvals[eig];
            r[i] = vec {M(i,0), M(i,1), M(i,2)};
        }
        v[eig][0] = (r[0]+r[1])[1]*r[2][2] - (r[0]+r[1])[2]*r[2][1];
        v[eig][1] = (r[0]+r[1])[2]*r[2][0] - (r[0]+r[1])[0]*r[2][2];
        v[eig][2] = (r[0]+r[1])[0]*r[2][1] - (r[0]+r[1])[1]*r[2][0];
    }
    return v;
}