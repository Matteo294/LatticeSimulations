#include "MMC.h"

MMC::MMC(class Model* s, class Lattice* l) : Simulator(s, l){;}

void MMC::runMC(int n, int Nskip, double thermalization){
    int accepted=0;
    int counter=0;
    double S, dS, oldS, avgS=0.;
    double r;
    S = s->evaluateAction();
    for(int i=0; i<n; i++){
        oldS = S;
        s->newConf();
        S = s->evaluateAction();
        dS = S-oldS;
        r = (double)rand()/RAND_MAX;
        if ((dS>0) && (exp(-dS) < r)){s->writeConfiguration();}
        else accepted++;
        if(n%Nskip == 0){avgS += S; counter++;}
    }
    cout << "Avg S: " << (double) avgS/counter << endl;
}
