#include "Phi4_2d.h"

Phi4_2d::Phi4_2d(class Lattice* latt, double m2, double g) : Model(latt){
    this->m2=m2;
    this->g=g;
    for(int nt=0; nt<lattice->Nt; nt++){
        for(int nx=0; nx<lattice->Nx; nx++){
            for(int nx=0; nx<lattice->Nx; nx++){
                for(int ny=0; ny<lattice->Ny; ny++){
                    for(int nz=0; nz<lattice->Nz; nz++){
                        phi[nt][nx][ny][nz] = gaussian(seed_gaussian);
                    }
                }
            }
        }
    }
}

double Phi4_2d::evaluateAction(){
    double sum = 0.0;
    int Nt = lattice->Nt;
    int Nx = lattice->Nx;
    int Ny = lattice->Ny;
    int Nz = lattice->Nz;
    for(int nt=0; nt<Nt; nt++){
        for(int nx=0; nx<Nx; nx++){
            for(int ny=0; ny<Ny; ny++){
                for(int nz=0; nz<Nz; nz++){
                    sum = sum - 0.5 * phi[nt][nx][ny][nz] * (
                        phi[PBCidx(nt+1, Nt)][nx][ny][nz] + phi[PBCidx(nt-1, Nt)][nx][ny][nz] + 
                        phi[nt][PBCidx(nx+1, Nx)][ny][nz] + phi[nt][PBCidx(nx-1, Nx)][ny][nz] + 
                        phi[nt][nx][PBCidx(ny+1, Ny)][nz] + phi[nt][nx][PBCidx(ny-1, Ny)][nz] +
                        phi[nt][nx][ny][PBCidx(nz+1, Nz)] + phi[nt][nx][ny][PBCidx(nz-1, Nz)] -
                        8*phi[nt][nx][ny][nz]
                        )
                        + 0.5*m2*phi[nt][nx][ny][nz]*phi[nt][nx][ny][nz] + (double) g/24. * pow(phi[nt][nx][ny][nz], 4);
                }
            }
        }
    }
    return sum;
}

void Phi4_2d::modelInfo(){
    cout << endl << "phi4 theory --> -phi*ddphi + m2*phi^2 + g*phi^4" << endl;
    cout << "Fields --> 1 scalar field phi" << endl;
    cout << "Couplings --> m2: " << m2 << " g: " << g << endl << endl;
}

double Phi4_2d::evaluateMDdrift(int nt, int nx, int ny, int nz){
    int Nt = lattice->Nt;
    int Nx = lattice->Nx;
    int Ny = lattice->Ny;
    int Nz = lattice->Nz;
    return 
        phi[PBCidx(nt+1, Nt)][nx][ny][nz] + phi[PBCidx(nt-1, Nt)][nx][ny][nz] + 
        phi[nt][PBCidx(nx+1, Nx)][ny][nz] + phi[nt][PBCidx(nx-1, Nx)][ny][nz] +
        phi[nt][nx][PBCidx(ny+1, Ny)][nz] + phi[nt][nx][PBCidx(ny-1, Ny)][nz] +
        phi[nt][nx][ny][PBCidx(nz+1, Nz)] + phi[nt][nx][ny][PBCidx(nz-1, Nz)] - 
        8*phi[nt][nx][ny][nz]
        - m2*phi[nt][nx][ny][nz] // phi^2 term
        - (double) g/6.*pow(phi[nt][nx][ny][nz], 3); // phi^4 term
}

vector<vector<vector<vector<double>>>> Phi4_2d::copyConfiguration(){
    return phi;
}

void Phi4_2d::writeConfiguration(vector<vector<vector<vector<double>>>> phi){
    for(int nt = 0; nt<lattice->Nt; nt++){
        for(int nx = 0; nx<lattice->Nx; nx++){
            for(int ny=0; ny<lattice->Ny; ny++){
                for(int nz=0; nz<lattice->Nz; nz++){
                    this->phi[nt][nx][ny][nz] = phi[nt][nx][ny][nz];
                }
            }
        }
    }
}