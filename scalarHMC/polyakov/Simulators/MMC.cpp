#include "MMC.h"

MMC::MMC(class Model* s, class Lattice* l) : Simulator(s, l){;}

void MMC::runMC(int NMC, int Nskip, double thermalization){
    srand(time(NULL));
    int accepted=0;
    int counter=0;
    double S, dS, oldS, avgS=0., avgS2=0.;
    double avgdS = 0.;
    double r;
    const int Ntherm = (int) ((double)thermalization*NMC);
    S = s->evaluateAction();
    for(int i=0; i<NMC; i++){
        for(int idx=0; idx<10; idx++){
            oldS = S;
            s->saveConfiguration();
            s->newConf(idx);
            S = s->evaluateAction();
            dS = S-oldS;
            avgdS += dS;
            r = (double)rand()/RAND_MAX;
            if ((dS>0) && (exp(-dS) < r)){s->writeSavedConfiguration(); S=oldS;}
            if((i%Nskip == 0) && (i>Ntherm)){avgS+=S; avgS2+=S*S; counter++;}
        }
    }
    avgS = (double) avgS/counter;
    avgS2 = (double) avgS2/counter;
    double stdDev = sqrt(avgS2 - avgS*avgS);
    cout << "Avg S: " <<  avgS << " StdDev: " << stdDev << endl;
    cout << (double)avgdS/counter << endl;
}
